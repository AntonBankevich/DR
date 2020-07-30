//
// Created by anton on 17.07.2020.
//
#define _GLIBCXX_PARALLEL
#include <iostream>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>
#include <common/dir_utils.hpp>
#include <queue>
#include <omp.h>
#include "common/cl_parser.hpp"
#include "rolling_hash.hpp"
#include "sparse_dbg.hpp"
#include<parallel/algorithm>
#include <sstream>
#include <unordered_set>
#include "hash_utils.hpp"
#include "common/omp_utils.hpp"
#include "common/logging.hpp"
#include "common/bloom_filter.hpp"
#include "common/output_utils.hpp"

typedef unsigned __int128 htype128;

size_t bit_count(unsigned char mask) {
    return ((mask >> 0u)& 1u) + ((mask >> 1u)& 1u) + ((mask >> 2u)& 1u) + ((mask >> 3u)& 1u) +
            ((mask >> 4u)& 1u) + ((mask >> 5u)& 1u) + ((mask >> 6u)& 1u) + ((mask >> 7u)& 1u);
}

std::vector<htype128> constructMinimizers(Time &time, const string &reads_file, size_t threads, const RollingHash<htype128> &hasher, const size_t w) {
    std::cout << "Reading reads" << std::endl;
    std::vector<std::vector<htype128>> prev;
    prev.resize(threads);
    const size_t buffer_size = 1000000000;
    std::cout << time.get() << "Extracting minimizers" << std::endl;
    size_t total_len = 0;
    size_t read_num = 0;
    std::vector<htype128> hash_list;
    io::SeqReader reader(reads_file);
    while(not reader.eof()) {
        size_t tlen = 0;
        std::cout << time.get() << "Starting new round" << std::endl;
        std::vector<Contig> reads;
        reads.reserve(1000000);
        std::vector<std::vector<htype128>> hashs;
        hashs.resize(threads);
#pragma omp parallel default(none) shared(hasher, prev, hashs, hash_list, reader, tlen, buffer_size, std::cout, time, reads, w)
        {
#pragma omp single
            {
                while (not reader.eof() and tlen < buffer_size) {
                    reads.push_back(reader.read());
                    tlen += reads.back().size();
                    if(reads.back().size() < hasher.k + w) {
                        continue;
                    }
                    size_t index = reads.size() - 1;
#pragma omp task default(none) shared(reads, hasher, hashs, std::cout, index, w)
                    {
                        const Contig &read = reads[index];
                        MinimizerCalculator<htype128> calc(read.seq, hasher, w);
                        std::vector<htype128> minimizers(calc.minimizerHashs());
                        std::vector<htype128> & thread_hashs = hashs[omp_get_thread_num()];
                        thread_hashs.insert(thread_hashs.end(), minimizers.begin(), minimizers.end());
                    }
                }
                for(std::vector<htype128> & tmp: prev) {
                    hash_list.insert(hash_list.end(), tmp.begin(), tmp.end());
                }
                std::cout << time.get() << tlen  << " nucleotides in " << reads.size() <<
                                " sequences were read from disk. Processing in progress  " << std::endl;
            }
        }
        std::swap(hashs, prev);
        total_len += tlen;
        read_num += reads.size();
        std::cout << time.get() << "Processing finished. Total nucleotides processed: " << total_len << std::endl;
        std::cout << time.get() << "Total reads processed: " << read_num << std::endl;
    }
    for(std::vector<htype128> & tmp: prev) {
        hash_list.insert(hash_list.end(), tmp.begin(), tmp.end());
    }
    std::cout << time.get() << "Finished read processing" << std::endl;
    std::cout << time.get() << hash_list.size() << " hashs collected. Starting sorting." << std::endl;
    //    TODO replace with parallel std::sort
    __gnu_parallel::sort(hash_list.begin(), hash_list.end());
    hash_list.erase(std::unique(hash_list.begin(), hash_list.end()), hash_list.end());
    std::cout << time.get() << "Finished sorting. Total distinct minimizers: " << hash_list.size() << std::endl;
    return hash_list;
}

std::vector<Sequence> extractDisjointigs(const Time & time, SparseDBG<htype128> &sdbg, size_t threads) {
    std::cout << time.get() << "Starting to extract disjointigs." << std::endl;
    std::vector<std::vector<Sequence>> res;
    res.resize(threads);
#pragma omp parallel default(none) shared(std::cout, time, sdbg, res)
    {
#pragma omp single
        {
            for(auto it = sdbg.v.begin(); it != sdbg.v.end(); ++it) {
                ExtensionsRecord<htype128> &rec = it->second;
                htype128 hash = it->first;
                if(!rec.isJunction())
                    continue;
#pragma omp task default(none) shared(sdbg, rec, res, hash)
                {
                    for(size_t i = 0; i < rec.getRightExtensions().size(); i++) {
                        res[omp_get_thread_num()].push_back(sdbg.walkForward(hash, rec, i));
                    }
                }
            }
        }
    }
    std::vector<Sequence> rres;
    for(const std::vector<Sequence> &r : res) {
        rres.insert(rres.end(), r.begin(), r.end());
    }
    size_t sz = 0;
    for (Sequence & seq: rres) {
        sz += seq.size();
//        std::cout << seq << std::endl;
    }
    std::cout << time.get() << "Finished extracting " << rres.size() << " disjointigs of total size " << sz << std::endl;
    return rres;
}

bool checkJunction(unsigned char mask) {
    unsigned char right = mask & 15u;
    unsigned char left = mask >> 4u;
    return bit_count(left) != 1 || bit_count(right) != 1;
}

template<class T>
size_t total_size(const std::vector<T> &data) {
    size_t res = 0;
    for(const T & val : data) {
        res += val.size();
    }
    return res;
}

std::vector<htype128> findJunctions(const Time time, const std::vector<Sequence>& disjointigs,
        const RollingHash<htype128> &hasher, size_t threads) {
    bloom_parameters parameters;
    parameters.projected_element_count = total_size(disjointigs);
    parameters.false_positive_probability = 0.0001;
    VERIFY(!!parameters);
    parameters.compute_optimal_parameters();
    BloomFilter filter(parameters);
    const RollingHash<htype128> ehasher = hasher.extensionHash();
    std::cout << time.get() << "Filling bloom filter with k+1-mers." << std::endl;
#pragma omp parallel for default(none) shared(filter, disjointigs, ehasher, hasher, std::cout)
    for(size_t i = 0; i < disjointigs.size(); i++) {
        const Sequence &seq = disjointigs[i];
        KWH<htype128> kmer(ehasher, seq, 0);
//        KWH<htype128> kmer1(hasher, seq, 0);
        while(true) {
            filter.insert(kmer.hash);
            if(!kmer.hasNext())
                break;
            kmer = kmer.next();
//            kmer1 = kmer1.next();
//            std::cout << kmer1.extendRight(seq[kmer1.pos + hasher.k]) << " "
//                      << kmer1.extendLeft(seq[kmer1.pos - 1]) << " " << kmer.hash << std::endl;
        }
    }
    std::cout << time.get() << "Finished filling bloom filter. Selecting junctions." << std::endl;
    ParallelRecordCollector<htype128> junctions(threads);
#pragma omp parallel for default(none) shared(filter, disjointigs, hasher, junctions, std::cout)
    for(size_t i = 0; i < disjointigs.size(); i++) {
        const Sequence &seq = disjointigs[i];
        KWH<htype128> kmer(hasher, seq, 0);
        while(true) {
            size_t cnt1 = 0;
            size_t cnt2 = 0;
            for(char c = 0; c < 4u; c++) {
                cnt1 += filter.contains(kmer.extendRight(c));
                cnt2 += filter.contains(kmer.extendLeft(c));
            }
            if(cnt1 != 1 || cnt2 != 1) {
                junctions.emplace_back(kmer.hash);
//                std::cout << kmer.seq.Subseq(kmer.pos, kmer.pos + hasher.k).str() << std::endl;
            }
            if(!kmer.hasNext())
                break;
            kmer = kmer.next();
        }
    }

    std::vector<htype128> res = junctions.collect();
    __gnu_parallel::sort(res.begin(), res.end());
    res.erase(std::unique(res.begin(), res.end()), res.end());
    std::cout << time.get() << "Collected " << res.size() << " junctions." << std::endl;
    return res;
}

std::vector<Sequence> constructDBG(const Time time, const std::vector<htype128> &vertices, const std::vector<Sequence> &disjointigs,
        const RollingHash<htype128> &hasher) {
    std::cout << time.get() << "Starting DBG construction." << std::endl;
    std::unordered_set<htype128> vset(vertices.begin(), vertices.end());
    std::unordered_map<htype128, std::vector<Sequence>> edges;
    for(auto & v: vertices) {
        edges[v] = std::vector<Sequence>(4);
    }
    std::vector<Sequence> res;
    for(const Sequence &seq: disjointigs) {
        KWH<htype128> kmer(hasher, seq, 0);
        std::vector<std::pair<htype128, size_t>> junctions;
        while(true) {
            if (vset.find(kmer.hash) != vset.end()) {
                junctions.emplace_back(kmer.hash, kmer.pos);
            }
            if(!kmer.hasNext())
                break;
            kmer = kmer.next();
        }
        for(size_t i = 0; i + 1 < junctions.size(); i++) {
            res.push_back(seq.Subseq(junctions[i].second, junctions[i + 1].second + hasher.k));
            edges[junctions[i].first][seq[junctions[i].second + hasher.k]] = res.back();
        }
    }
    size_t cnt = 0;
    for(auto & v: vertices)
        for(Sequence &seq: edges[v])
            if (!seq.empty())
                cnt += 1;
    std::cout << time.get() << "Constructed " << res.size() << " edge sequences with " << cnt << " unique."  << std::endl;
    return res;
}

int main(int argc, char **argv) {
    Time time;
    CLParser parser({"reads=", "unique=none", "output-dir=", "threads=8", "k-mer-size=7000", "window=3000", "base=239", "debug"},
            {"o=output-dir", "t=threads", "k=k-mer-size","w=window"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    const RollingHash<htype128> hasher(std::stoi(parser.getValue("k-mer-size")), std::stoi(parser.getValue("base")));
    const size_t w = std::stoi(parser.getValue("window"));
    std::string reads_file = parser.getValue("reads");
    size_t threads = std::stoi(parser.getValue("threads"));
    omp_set_num_threads(threads);

    std::vector<htype128> hash_list;
    if (parser.getValue("unique") != "none") {
        std::ifstream is;
        is.open(parser.getValue("unique"));
        readHashs(is, hash_list);
        is.close();
    } else {
        hash_list = constructMinimizers(time, reads_file, threads, hasher, w);
        std::ofstream os;
        os.open(std::string(dir.c_str()) + "/unique.save");
        writeHashs(os, hash_list);
        os.close();
    }

    SparseDBG<htype128> sdbg = constructSparseDBGFromReads(time, reads_file, threads, hasher, std::move(hash_list), w);
    sdbg.printStats();
//    sdbg.print();

    tieTips(time, sdbg, w, threads);
    sdbg.printStats();
//    sdbg.print();

    std::vector<Sequence> disjointigs = extractDisjointigs(time, sdbg, threads);
    std::vector<htype128> vertices = findJunctions(time, disjointigs, hasher, threads);
    std::vector<Sequence> dbg = constructDBG(time, vertices, disjointigs, hasher);

    return 0;
}