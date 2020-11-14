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

struct AlignmentStats {
    size_t no_align = 0;
    size_t multi_align = 0;
    size_t single_align = 0;
    size_t full_align = 0;
    std::vector<size_t> readErrors = std::vector<size_t>(100);
    std::vector<size_t> refErrors = std::vector<size_t>(4000000);
    std::vector<size_t> cov = std::vector<size_t>(4000000);

    AlignmentStats operator+(const AlignmentStats &other) const {
        AlignmentStats res;
        res.no_align = no_align + other.no_align;
        res.multi_align = multi_align + other.multi_align;
        res.single_align = single_align + other.single_align;
        res.full_align = full_align + other.full_align;
        for(size_t i = 0; i < readErrors.size(); i++) {
            res.readErrors[i] = readErrors[i] + other.readErrors[i];
        }
        for(size_t i = 0; i < refErrors.size(); i++) {
            res.refErrors[i] = refErrors[i] + other.refErrors[i];
        }
        return std::move(res);
    }
};

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

void FillStat(AlignmentStats &stat, PositionalAlignment<Contig, Contig> &best) {
    size_t diff = 0;
    size_t prev = 0;
    for (size_t i = 0; i + 1 < best.positions_from.size(); i++) {
        if (best.positions_from[i + 1] > best.positions_from[i] + 1 ||
            best.positions_to[i + 1] > best.positions_to[i] + 1) {
            stat.refErrors[best.positions_to[i] / 1000] += 1;
            diff += 1;
        }
    }
    for (size_t i = best.seg_to.left / 1000; i < best.seg_to.right / 1000; i++) {
        stat.cov[i] += 1;
    }
    stat.readErrors[std::min(diff, stat.readErrors.size() - 1)] += 1;
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
    CLParser parser({"reference=", "output-dir=", "threads=8"}, {"reads"},
                    {"o=output-dir", "t=threads"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    size_t threads = std::stoi(parser.getValue("threads"));
    omp_set_num_threads(threads);
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "analysis");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile());
    logger.info() << "Reading genome" << std::endl;
    StringContig::needs_compressing = true;
    std::vector<Contig> ref = io::SeqReader(parser.getValue("reference")).readAllContigs();
    SequenceBuilder sb;
    sb.appendAll(ref.begin(), ref.end());
    Contig concatRef(sb.BuildSequence(), "concat");
    ref.clear();
    ref.push_back(concatRef);
    ref.push_back(concatRef.RC());
    logger.info() << "Finished reading genome" << std::endl;
    logger.info() << "Building minimap index of reference" << std::endl;
    RawAligner<Contig> aligner(ref, threads, "ava-hifi");
    logger.info() << "Finished building minimap index" << std::endl;
    logger.info() << "Aligning and analysing reads" << std::endl;
    io::Library libReads = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::Library libCorrected = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("corrected"));
    io::SeqReader reads(libReads);
    ParallelRecordCollector<Segment<Contig>> segs(threads);

    std::function<void(size_t, const Contig &, const std::vector<CigarAlignment<Contig, Contig>> &)> collect_task=
            [&segs] (size_t ind, const Contig &read, const std::vector<CigarAlignment<Contig, Contig>> &cigar_als) {
                std::vector<PositionalAlignment<Contig, Contig>> als = selectAlignments(cigar_als);
                size_t best = chooseBest(als, read);
                if(best != size_t(-1)) {
                    segs.emplace_back(cigar_als[best].seg_to);
                    VERIFY(als[best].pi() > 0.5);
                }
            };
    alignment_recipes::AlignAndProcess(reads.begin(), reads.end(), aligner, collect_task, logger, threads);

    logger.info() << "Collecting events" << std::endl;
    std::vector<std::pair<size_t, int>> events;
    for(Segment<Contig> & seg : segs) {
        if(seg.contig() == concatRef){
            events.emplace_back(seg.left, 1);
            events.emplace_back(seg.right, -1);
        } else {
            events.emplace_back(concatRef.size() - seg.right, 1);
            events.emplace_back(concatRef.size() - seg.left, -1);
        }
    }
    logger.info() << "Sorting events" << std::endl;
    std::sort(events.begin(), events.end(), [] (const std::pair<size_t, int> &e1, const std::pair<size_t, int> &e2)
              {
                    return e1.first < e2.first || (e1.first == e2.second && e1.second > e2.first);
              }
    );
    logger.info() << "Printing result" << std::endl;
    size_t start = 0;
    int cov = 0;
    int max_cov = 0;
    std::vector<std::pair<double, Segment<Contig>>> bad;
    for(size_t i = 0; i < events.size(); i++) {
        cov +=events[i].second;
        max_cov = std::max(max_cov, cov);
        if(cov == 0) {
            if(max_cov > 3) {
                double ratio = ( i - start + 1.0) / std::max(events[i].first - events[start].first - 8000, size_t(1));
                bad.emplace_back(-ratio, concatRef.segment(events[start].first, events[i].first));
                logger << "New seg " << events[start].first << " " << events[i].first - events[start].first << " " << (i - start + 1) / 2<< std::endl;
                for(size_t j = start; j <= i; j++) {
                    logger << events[j].first - events[start].first << " " << events[j].second << std::endl;
                }
            } else {
                logger << "Low covered seg " << events[start].first << " " << events[i].first - events[start].first << std::endl;
            }
            start = i + 1;
            max_cov = 0;
        }
    }
    std::sort(bad.begin(), bad.end());
    std::ofstream os;
    os.open(dir / "segments.fasta");
    for(std::pair<double, Segment<Contig>> & pair : bad) {
        os << ">" << pair.second.left << "-" << pair.second.size() << "-" << -pair.first << std::endl;
        os << pair.second.seq() << std::endl;
    }
    os.close();
    return 0;
}