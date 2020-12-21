//
// Created by Andrey Bzikadze on 12/17/20.
//

#include <iostream>
#include <common/cl_parser.hpp>
#include <sequences/seqio.hpp>
#include "../dbg/rolling_hash.hpp"
#include "../dbg/hash_utils.hpp"

using UMapSeq = std::unordered_map<std::string, Sequence>;

template<typename T>
using UMapSeq2KWH = std::unordered_map<std::string, KWH<T>>;

using PosUMap = std::unordered_map<std::string, std::vector<std::pair<std::string, int>>>;

UMapSeq read_fasta(const std::vector<std::string>& path,
                   bool verbose=false) {
    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(path);
    io::SeqReader seq_reader(lib);

    UMapSeq seqs;
    for(StringContig tmp : seq_reader) {
        if (verbose) {
            std::cout << tmp.id << " " << tmp.size() << std::endl;
        }
        seqs.emplace(tmp.id, Sequence(tmp.seq));
    }
    return seqs;
}

size_t get_min_length(const UMapSeq & seqs) {
    size_t min_len = std::numeric_limits<size_t>::max();
    for (const auto& pair: seqs) {
        auto& seq = pair.second;
        min_len = std::min(min_len, seq.size());
    }
    return min_len;
}

template<typename T>
UMapSeq2KWH<T> get_KWH(const UMapSeq& seqs, const RollingHash<T>& hasher) {
    UMapSeq2KWH<T> hashes;
    for (const auto& pair: seqs) {
        auto & s_id = pair.first;
        auto & seq = pair.second;
        hashes.emplace(s_id, KWH<T>(hasher, seq, 0));
    }
    return hashes;
}

template<typename T>
std::pair<PosUMap, PosUMap>
find_matches(const UMapSeq2KWH<T>& queries_kwh, const UMapSeq2KWH<T>& targets_kwh) {
    std::unordered_map<T, std::vector<std::string>> hash2seqs;
    for (const auto& pair: queries_kwh) {
        const auto& s_id = pair.first;
        const auto& kwh = pair.second;
        T hash = kwh.hash();
        hash2seqs[hash].emplace_back(s_id);
    }

    PosUMap pos_f, pos_r;
    for (const auto& pair: queries_kwh) {
        const auto& q_id = pair.first;
        pos_f[q_id];
        pos_r[q_id];
    }

    for (const auto& pair: targets_kwh) {
        const auto& t_id = pair.first;
        std::cout << t_id << std::endl;
        const auto& target_kwh = pair.second;

        KWH<T> kwh(target_kwh);
        while(true) {
            if (kwh.pos % 10000000 == 0) {
                std::cout << kwh.pos << std::endl;
            }
            const T hash = kwh.hash();
            const auto search = hash2seqs.find(hash);
            if (search != hash2seqs.end()) {
                const auto& q_ids = search->second;
                for (const auto& q_id: q_ids) {
                    const auto& query = queries_kwh.at(q_id).getSeq();
                    const auto target = kwh.getSeq();
                    if (query == target) {
                        pos_f[q_id].emplace_back(std::make_pair(t_id, kwh.pos));
                    }
                    if (!query == target) {
                        pos_r[q_id].emplace_back(std::make_pair(t_id, kwh.pos));
                    }
                }
            }
            if (!kwh.hasNext()) {
                break;
            }
            kwh = kwh.next();
        }
    }
    return { pos_f, pos_r };
}

void output_matches(const PosUMap& pos_f, const PosUMap& pos_r, const std::string& path) {
    const std::experimental::filesystem::path outfn(path);
    std::ofstream os;
    os.open(outfn);
    os << "query_id\tmult_forward\tmult_reverse\n";
    for (const auto& pair : pos_f) {
        const auto q_id = pair.first;
        const auto& q_pos_f = pair.second;
        size_t cnt_q_pos_f = q_pos_f.size();
        const auto& q_pos_r = pos_r.at(q_id);
        size_t cnt_q_pos_r = q_pos_r.size();
        os << q_id << "\t" << cnt_q_pos_f << "\t" << cnt_q_pos_r << std::endl;
    }
}

int main(int argc, char** argv) {
    CLParser parser({"output=", "reference=", "query=", "base=239"}, {},
                    {"o=output", "r=reference", "q=query"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }

    auto ref = read_fasta(parser.getListValue("reference"), true);
    auto queries = read_fasta(parser.getListValue("query"));

    size_t q_min_len = get_min_length(queries);
    std::cout << "Queries minimum length " << q_min_len << std::endl;

    const int base = std::stoi(parser.getValue("base"));

    RollingHash<htype128> hasher(q_min_len, base); // this object has to live til the end
    const auto queries_kwh = get_KWH<htype128>(queries, hasher);
    const auto ref_kwt = get_KWH<htype128>(ref, hasher);

    const auto pair_pos = find_matches(queries_kwh, ref_kwt);
    const auto& pos_f = pair_pos.first;
    const auto& pos_r = pair_pos.second;

    output_matches(pos_f, pos_r, parser.getValue("output"));
}