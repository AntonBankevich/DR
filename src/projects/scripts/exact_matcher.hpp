#pragma once
//
// Created by Andrey Bzikadze on 12/21/20.
//

#include <iostream>
#include <sequences/seqio.hpp>
#include "../dbg/rolling_hash.hpp"
#include "../dbg/hash_utils.hpp"

namespace scripts {
namespace exact_matcher {

enum class Strand {
    forward, reverse
};

struct QueryPositionInTarget {
    const Contig &target;
    const Contig &query;
    const size_t pos;
    const Strand strand;

    QueryPositionInTarget(const Contig &target_, const Contig &query_, const size_t pos_,
                          const Strand strand_) :
            target{target_}, query{query_}, pos{pos_}, strand{strand_} {
        VERIFY(pos_ + query.size() <= target.size())
    }
};

using PositionsVector = std::vector<scripts::exact_matcher::QueryPositionInTarget>;

namespace details {

    size_t get_min_length(const std::vector<Contig> &contigs) {
        size_t min_len = std::numeric_limits<size_t>::max();
        for (const Contig &contig: contigs) {
            min_len = std::min(min_len, contig.size());
        }
        return min_len;
    }

    template<typename htype>
    struct NamedKWH {
        const Contig &contig;
        KWH<htype> kwh;

        NamedKWH(const Contig &contig_, KWH<htype> kwh_) : contig{contig_}, kwh{kwh_} {}
    };

    template<typename htype>
    std::vector<NamedKWH<htype>>
    get_named_KWH(const std::vector<Contig> &contigs, const RollingHash<htype> &hasher) {
        std::vector<NamedKWH<htype>> named_kwh;
        for (const Contig &contig: contigs) {
            named_kwh.emplace_back(contig, KWH<htype>(hasher, contig.seq, 0));
        }
        return named_kwh;
    }


    template<typename htype>
    PositionsVector find_matches(const std::vector<NamedKWH<htype>> &queries_kwh,
                                 const std::vector<NamedKWH<htype>> &targets_kwh) {
        PositionsVector positions;
        std::unordered_map<htype, const Contig *> hash2seqs;
        for (const NamedKWH<htype> &query_kwh: queries_kwh) {
            const Contig *ctg_pointer = &query_kwh.contig;
            const KWH<htype> &kwh = query_kwh.kwh;
            htype hash = kwh.hash();
            hash2seqs.emplace(hash, ctg_pointer);
        }

        for (const NamedKWH<htype> &target_kwh : targets_kwh) {
            KWH<htype> rolling_target_kwh(target_kwh.kwh);
            while (true) {
                const htype hash = rolling_target_kwh.hash();
                const auto search = hash2seqs.find(hash);
                if (search != hash2seqs.end()) {
                    for (const NamedKWH<htype> &query_kwh : queries_kwh) {
                        const Contig &query = query_kwh.contig;
                        size_t pos = rolling_target_kwh.pos;
                        if (pos + query.size() > target_kwh.contig.size()) {
                            continue;
                        }
                        const Sequence target_substr = target_kwh.contig.seq.Subseq(pos, pos + query.size());
                        if (query.seq == target_substr) {
                            positions.emplace_back(target_kwh.contig,
                                                   query_kwh.contig,
                                                   pos,
                                                   Strand::forward);
                        }
                        if (query.seq == !target_substr) {
                            positions.emplace_back(target_kwh.contig,
                                                   query_kwh.contig,
                                                   pos,
                                                   Strand::reverse);
                        }
                    }
                }
                if (!rolling_target_kwh.hasNext()) {
                    break;
                }
                rolling_target_kwh = rolling_target_kwh.next();
            }
        }
        return positions;
    }

} // End namespace details

template<typename htype>
PositionsVector exact_match(const std::string &queries_path,
                            const std::string &targets_path,
                            const size_t base = 239) {
    using namespace details;
    io::SeqReader queriesReader(queries_path);
    io::SeqReader targetsReader(targets_path);

    const std::vector<Contig> queries = queriesReader.readAllContigs();
    const std::vector<Contig> targets = targetsReader.readAllContigs();

    if (queries.empty() or targets.empty()) {
        return PositionsVector();
    }

    const size_t query_min_len = details::get_min_length(queries);

    const RollingHash<htype> hasher(query_min_len, base); // this object has to live till the end

    const std::vector<NamedKWH<htype>> queries_kwh = get_named_KWH<htype>(queries, hasher);
    const std::vector<NamedKWH<htype>> targets_kwh = get_named_KWH<htype>(targets, hasher);

    PositionsVector positions = find_matches(queries_kwh, targets_kwh);

    return positions;
}

void output_matches(const std::vector<QueryPositionInTarget> &positions, const std::string &out_path) {
    using namespace details;
    const std::experimental::filesystem::path out_fn(out_path);
    std::ofstream os;
    os.open(out_fn);
    os << "target_id\tquery_id\tposition\tstrand\n";
    for (const QueryPositionInTarget &query_pos_in_target : positions) {
        const std::string query_id = query_pos_in_target.query.id;
        const std::string target_id = query_pos_in_target.target.id;
        const size_t position = query_pos_in_target.pos;
        const Strand strand = query_pos_in_target.strand;
        const std::string strand_seq = strand == Strand::forward ? "+" : "-";
        os << target_id << '\t' << query_id << '\t' << position << '\t' << strand_seq << std::endl;
    }
}

} // End namespace exact_matcher
} // End namespace scripts