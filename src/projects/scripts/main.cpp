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
    CLParser parser({"reference=", "output-dir=", "threads=8"}, {"reads", "corrected"},
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
    std::vector<Contig> ref = io::SeqReader(parser.getValue("reference")).readAllCompressedContigs();
    SequenceBuilder sb;
    sb.appendAll(ref.begin(), ref.end());
    Contig concatRef(sb.BuildSequence(), "concat");
    ref.clear();
    ref.push_back(concatRef);
    ref.push_back(concatRef.RC());
    logger.info() << "Finished reading genome" << std::endl;
    logger.info() << "Building minimap index of reference" << std::endl;
    RawAligner<Contig> aligner(ref, threads, "ava-hifi");
    logger.info() << "Aligner parts " << aligner.index.size() << std::endl;
    logger.info() << "Finished building minimap index" << std::endl;
    logger.info() << "Aligning and analysing reads" << std::endl;
    io::Library libReads = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::Library libCorrected = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("corrected"));
    io::SeqReader reads(libReads);
    io::SeqReader corrected(libCorrected);
    std::vector<AlignmentStats> stats(threads);
    std::vector<AlignmentStats> stats_corrected(threads);
    ParallelRecordCollector<Contig> strange_reads(threads);
    ParallelRecordCollector<Contig> bad_reads(threads);
    ParallelRecordCollector<Contig> nonperfect_reads(threads);
    ParallelRecordCollector<Contig> unaligned_reads(threads);
    ParallelRecordCollector<Contig> perfect_reads(threads);
//    std::vector<Contig> contigs;
//    for(StringContig contig : corrected) {
//        contigs.push_back(contig.makeCompressedContig());
//    }
//    logger.info() << "Batch test" << std::endl;
//    aligner.align1(contigs, threads);
//    logger.info() << "Sequential test" << std::endl;
//#pragma omp parallel for default(none) shared(aligner, contigs)
//    for(size_t i = 0; i < contigs.size(); i++) {
//        aligner.align(contigs[i]);
//    }
//    logger.info() << "Done" << std::endl;
//    corrected.reset();
    std::vector<double> pis;
    std::vector<Contig> initial;
    pis.resize(10000000);
    initial.resize(10000000);

    std::function<void(size_t, const Contig &, const std::vector<CigarAlignment<Contig, Contig>> &)> collect_task=
        [&pis, &initial, &unaligned_reads] (size_t ind, const Contig &read, const std::vector<CigarAlignment<Contig, Contig>> &cigar_als) {
            std::vector<PositionalAlignment<Contig, Contig>> als = selectAlignments(cigar_als);
            initial[ind] = read;
            size_t best = chooseBest(als, read);
            if(best == size_t(-1)) {
                pis[ind] = 0;
            } else {
                pis[ind] = als[best].pi();
                VERIFY(als[best].pi() > 0.5);
            }
//            if(best == size_t(-1)) {
//#pragma omp critical
//                {
//                    std::cout << "Unaligned " << read.id << std::endl;
//                }
//            }
        };
    std::function<void(size_t, const Contig &, const std::vector<CigarAlignment<Contig, Contig>> &)> compare_task =
            [&pis, &initial, &perfect_reads, &nonperfect_reads, &bad_reads, &strange_reads, &unaligned_reads]
                        (size_t ind, const Contig &corrected, const std::vector<CigarAlignment<Contig, Contig>> &cigar_als){
                std::vector<PositionalAlignment<Contig, Contig>> als = selectAlignments(cigar_als);
                size_t best = chooseBest(als, corrected);
                std::vector<std::string> tmp = split(corrected.id);
                size_t score = 0;
                if(tmp.size() >= 4) {
                    score = std::stoll(tmp[3]);
                }
                if(als.empty()) {
                    unaligned_reads.emplace_back(initial[ind]);
                } else if (als[best].isPerfect()) {
                    perfect_reads.emplace_back(corrected);
                } else {
                    nonperfect_reads.emplace_back(corrected);
                    if(score > 25000) {
                        bad_reads.emplace_back(initial[ind]);
                    } else if (score > 0 && score != size_t(-1)) {
                        strange_reads.emplace_back(initial[ind]);
                    }
                }
                if(!als.empty()) {
                    if (!als[best].isPerfect() || pis[ind] < 1) {
#pragma omp critical
                        {
                            std::cout << size_t((1 - pis[ind]) * 100000) * 0.001 << " " <<
                                      size_t((1 - als[best].pi()) * 100000) * 0.001 << " " <<
                                      corrected.id.substr(initial[ind].id.size()) << std::endl;
                        }
                    }
                } else {
#pragma omp critical
                    {
                        std::cout << "Fail " << als.size() << " " << cigar_als.size() << " " << pis[ind] << std::endl;
                    }
                }

            };

    alignment_recipes::AlignAndProcess(reads.begin(), reads.end(), aligner, collect_task, logger, threads);
    alignment_recipes::AlignAndProcess(corrected.begin(), corrected.end(), aligner, compare_task, logger, threads);

    logger.info() << "Not aligned " << unaligned_reads.size() << std::endl;
    logger.info() << "Perfect " << perfect_reads.size() << std::endl;
    logger.info() << "Nonperfect " << nonperfect_reads.size() << std::endl;
    logger.info() << "Strange " << strange_reads.size() << std::endl;
    logger.info() << "Bad " << bad_reads.size() << std::endl;
    PrintContigs(unaligned_reads, dir / "unaligned.fasta");
    PrintContigs(perfect_reads, dir / "perfect.fasta");
    PrintContigs(nonperfect_reads, dir / "nonperfect.fasta");
    PrintContigs(strange_reads, dir / "strange.fasta");
    PrintContigs(bad_reads, dir / "bad.fasta");
    logger.info() << "Printed contigs to files" << std::endl;
//    AlignmentStats res;
//    for(const AlignmentStats & stat : stats) {
//        res = res + stat;
//    }
//    logger.info() << "No read alignment " << res.no_align << std::endl;
//    logger.info() << "Single read alignment " << res.single_align << std::endl;
//    logger.info() << "Single read alignment is end-to-end " << res.full_align << std::endl;
//    logger.info() << "Multiple read alignment " << res.multi_align << std::endl;
//    for(size_t i = 0; i < res.readErrors.size(); i++) {
//        logger.info() << i << " " << res.readErrors[i] << std::endl;
//    }
//    std::ofstream os;
//    os.open(dir / "ref_errors.info");
//    for(size_t i = 0; i < res.cov.size(); i++) {
//        os << res.refErrors[i] << " " << res.cov[i] << std::endl;
//    }
//    os.close();
    return 0;
}