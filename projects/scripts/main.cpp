//
// Created by anton on 8/7/20.
//

#include <alignment/aligner.hpp>
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <vector>
#include <sequences/seqio.hpp>

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

int main(int argc, char **argv) {
    CLParser parser({"reads=", "reference=", "output-dir=", "threads=8"},
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
    logger << "Reading genome" << std::endl;
    std::vector<Contig> ref = io::SeqReader::CompressingReader(parser.getValue("reference")).readAll();
    SequenceBuilder sb;
    sb.appendAll(ref.begin(), ref.end());
    Contig concatRef(sb.BuildSequence(), "concat");
    ref.clear();
    ref.push_back(concatRef);
    logger << "Finished reading genome" << std::endl;
    logger << "Building minimap index of reference" << std::endl;
    RawAligner<Contig> aligner(ref, threads);
    logger << "Finished building minimap index" << std::endl;
    logger << "Aligning and analysing reads" << std::endl;
    io::SeqReader reader(io::SeqReader::CompressingReader(parser.getValue("reads")));
    std::vector<AlignmentStats> stats(threads);
#pragma omp parallel default(none) shared(reader, aligner, stats, logger)
    {
#pragma omp single
        {
            size_t cnt = 0;
            while(!reader.eof()) {
                cnt += 1;
                if (cnt % 100000 == 0)
                    logger << "Starting analysis of read number " << cnt << std::endl;
                Contig read = reader.read();
#pragma omp task default(none) shared(read, aligner, stats)
                {
                    std::vector<CigarAlignment<Contig, Contig>> cigar_als = aligner.align(read);
                    std::vector<PositionalAlignment<Contig, Contig>> als =
                            oneline::initialize<PositionalAlignment<Contig, Contig>>(cigar_als.begin(), cigar_als.end());
                    size_t min_al_size = std::max(std::max<size_t>(read.size(), 1000u) - 1000, read.size() * 3 / 4);
                    std::function<bool(PositionalAlignment<Contig, Contig>&)> f = [min_al_size](PositionalAlignment <Contig, Contig>& al) {
                        return al.seg_from.size() > min_al_size;
                    };
                    als = oneline::filter(als.begin(), als.end(), f);
                    AlignmentStats &stat = stats[omp_get_thread_num()];
                    if(als.empty()) {
                        stat.no_align += 1;
                    } else {
                        if (als.size() == 1) {
                            stat.single_align += 1;
                            if (als[0].seg_from.size() >= read.size() - 50) {
                                stat.full_align += 1;
                            }
                        } else {
                            stat.multi_align += 1;
                        }
                        PositionalAlignment<Contig, Contig>* best = &als[0];
                        for(PositionalAlignment<Contig, Contig>& al : als) {
                            if(best->seg_from.size() < read.size() - 50 && al.seg_from.size() >= read.size() - 50) {
                                best = &al;
                            }
                            if(best->seg_from.size() >= read.size() - 50 && al.seg_from.size() < read.size() - 50) {
                                continue;
                            }
                            if(best->pi() < al.pi()) {
                                best = &al;
                            }
                        }
                        size_t diff;
                        size_t prev = 0;
                        for(size_t i = 0; i + 1 < best->positions_from.size(); i++) {
                            if (best->positions_from[i + 1] > best->positions_from[i] + 1 || best->positions_to[i + 1] > best->positions_to[i] + 1){
                                stat.refErrors[best->positions_to[i] / 1000] += 1;
                                diff += 1;
                            }
                        }
                        for(size_t i = best->seg_to.left / 1000; i < best->seg_to.right / 1000; i++) {
                            stat.cov[i] += 1;
                        }
                        stat.readErrors[std::min(diff, stat.readErrors.size() - 1)];
                    }
                }
            }
        }
    }
    AlignmentStats res;
    for(const AlignmentStats & stat : stats) {
        res = res + stat;
    }
    logger << "No read alignment " << res.no_align << std::endl;
    logger << "Single read alignment " << res.single_align << std::endl;
    logger << "Single read alignment is end-to-end " << res.full_align << std::endl;
    logger << "Multiple read alignment " << res.multi_align << std::endl;
    for(size_t i = 0; i < res.readErrors.size(); i++) {
        logger << i << " " << res.readErrors[i] << std::endl;
    }
    std::ofstream os;
    os.open(dir / "ref_errors.info");
    for(size_t i = 0; i < res.cov.size(); i++) {
        os << res.refErrors[i] << " " << res.cov[i] << std::endl;
    }
    os.close();
    return 0;
}