//
// Created by anton on 8/7/20.
//

#include <alignment/aligner.hpp>
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <vector>

size_t stoi1(std::string s) {
    size_t res = 0;
    for(char c : s) {
        if(c < '0' || c > '9') {
            abort();
        }
        res = res * 10 + c - '0';
    }
    return res;
}

std::vector<PositionalAlignment<Contig, Contig>> selectAlignments(const std::vector<CigarAlignment<Contig, Contig>> &cigar_als) {
    std::vector<PositionalAlignment<Contig, Contig>> als =
            oneline::initialize<PositionalAlignment<Contig, Contig>>(cigar_als.begin(),
                                                                     cigar_als.end());
    if(als.empty()) {
        return als;
    }
    Contig read = cigar_als[0].seg_from.contig();
    size_t min_al_size = std::max(std::max<size_t>(read.size(), 1000u) - 1000, read.size() * 3 / 4);
    std::function<bool(PositionalAlignment<Contig, Contig> &)> f = [min_al_size](
            PositionalAlignment<Contig, Contig> &al) {
        return al.seg_from.size() > min_al_size;
    };
    als = oneline::filter(als.begin(), als.end(), f);
    return als;
}

size_t chooseBest(std::vector<PositionalAlignment<Contig, Contig>> &als, const Contig & read) {
    if (als.empty()) {
        return size_t(-1);
    } else {
        size_t res = 0;
        PositionalAlignment<Contig, Contig> *best = &als[0];
        size_t cnt = 0;
        for (PositionalAlignment<Contig, Contig> &al : als) {
            if (best->seg_from.size() < read.size() - 50 &&
                al.seg_from.size() >= read.size() - 50) {
                best = &al;
                res = cnt;
            }
            if (best->seg_from.size() >= read.size() - 50 &&
                al.seg_from.size() < read.size() - 50) {
                continue;
            }
            if (best->pi() < al.pi()) {
                best = &al;
                res = cnt;
            }
            cnt += 1;
        }
        return res;
    }
}

void PrintContigs(ParallelRecordCollector<Contig> &collector, const std::experimental::filesystem::path &path) {
    std::ofstream out;
    out.open(path);
    for(Contig &contig : collector) {
        out << ">" << contig.id << std::endl << contig.seq << std::endl;
    }
    out.close();
}

int main(int argc, char **argv) {
    CLParser parser({"reference=", "downsample=none", "threads=8"}, {"reads"},
                    {"t=threads"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    size_t threads = std::stoi(parser.getValue("threads"));
    omp_set_num_threads(threads);
    logging::Logger logger;
    logger << "Reading genome" << std::endl;
    std::vector<Contig> ref = io::SeqReader(parser.getValue("reference")).readAllCompressedContigs();
    SequenceBuilder sb;
    sb.appendAll(ref.begin(), ref.end());
    Contig concatRef(sb.BuildSequence(), "concat");
    ref.clear();
    ref.push_back(concatRef);
    ref.push_back(concatRef.RC());
    logger << "Finished reading genome" << std::endl;
    logger << "Building minimap index of reference" << std::endl;
    RawAligner<Contig> aligner(ref, threads, "ava-hifi");
    logger << "Aligner parts " << aligner.index.size() << std::endl;
    logger << "Finished building minimap index" << std::endl;
    logger << "Aligning and analysing reads" << std::endl;
    io::Library libReads = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::SeqReader reads(libReads);
    ParallelCounter perfect_reads(threads);
    ParallelCounter all_reads(threads);
    ParallelCounter len(threads);
    ParallelCounter errors(threads);
    ParallelRecordCollector<size_t> hist(threads);


    std::function<void(size_t, const Contig &, const std::vector<CigarAlignment<Contig, Contig>> &)> compare_task =
            [&perfect_reads, &all_reads, &len, &errors, &hist]
                    (size_t ind, const Contig &corrected, const std::vector<CigarAlignment<Contig, Contig>> &cigar_als){
                std::vector<PositionalAlignment<Contig, Contig>> als = selectAlignments(cigar_als);
                size_t best = chooseBest(als, corrected);
                if(!als.empty()){
                    if(als[best].isPerfect()) {
                        ++perfect_reads;
                        hist.emplace_back(0);
                    } else {
                        errors += als[best].diff_number();
                        hist.emplace_back(als[best].diff_number());
                    }
                }
                ++all_reads;
                len += corrected.size();
            };
    if(parser.getValue("downsample") == "none")
        alignment_recipes::AlignAndProcess(reads.begin(), reads.end(), aligner, compare_task, logger, threads);
    else {
        size_t rnum = std::stoull(parser.getValue("downsample"));
        std::vector<StringContig> readSeqs;
        for(StringContig read : reads) {
            readSeqs.emplace_back(std::move(read));
            if (readSeqs.size() == rnum)
                break;
        }
        alignment_recipes::AlignAndProcess(readSeqs.begin(), readSeqs.end(), aligner, compare_task, logger, threads);
    }
    logger << "Perfect " << perfect_reads.get() << std::endl;
    logger << "Errors " << errors.get() << std::endl;
    logger << "Errorrate " << double(errors.get()) / len.get() << std::endl;
    logger << "All " << all_reads.get() << std::endl;
    std::vector<size_t> h(30);
    for (size_t e : hist) {
        h[std::min(h.size() - 1, e)] += 1;
    }
    for(size_t i = 0; i < h.size(); i++) {
        std::cout << i << ": " << h[i] << std::endl;
    }
    return 0;
}