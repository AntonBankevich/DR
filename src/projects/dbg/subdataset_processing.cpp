#include "subdataset_processing.hpp"

#include <utility>

std::vector<RepeatResolver::Subdataset> RepeatResolver::SplitDataset(const std::function<bool(const Edge &)> &is_unique) {
    size_t k = dbg.hasher().getK();
    std::vector<Component> comps = ConditionSplitter(is_unique).splitGraph(dbg);
    std::vector<std::ofstream *> os;
    std::vector<std::ofstream *> alignments;
    std::vector<Subdataset> result;
    for(size_t i = 0; i < comps.size(); i++) {
        os.emplace_back(new std::ofstream());
        alignments.emplace_back(new std::ofstream());
        result.emplace_back(comps[i], dir / std::to_string(i));
        ensure_dir_existance(result.back().dir);
        os.back()->open(result.back().dir / "corrected.fasta");
        alignments.back()->open(result.back().dir / "alignments.txt");
        std::ofstream gos;
        gos.open(result.back().dir / "graph.fasta");
        printFasta(gos, result.back().component);
        gos.close();
        std::ofstream log;
        log.open(result.back().dir / "dbg.log");
        log << "-k " << k << std::endl;
        log.close();
        std::ofstream dot;
        dot.open(result.back().dir / "graph.dot");
        printDot(dot, result.back().component, readStorage.labeler());
        dot.close();
    }
    for(AlignedRead &read : readStorage) {
        if(!read.valid() || read.path.size() == 1)
            continue;
        GraphAlignment al = read.path.getAlignment();
        Sequence seq = al.Seq();
        std::stringstream ss;
        ss  << read.id << " " << al.start().hash() << int(al.start().isCanonical())
            << " " << read.path.cpath().str() << "\n";
        CompactPath rc(al.RC());
        ss  << "-" << read.id << " " << rc.start().hash() << int(rc.start().isCanonical())
            << " " << rc.cpath().str() << "\n";
        std::string alignment_record = ss.str();
        for(size_t j = 0; j < result.size(); j++) {
            for(size_t i = 1; i < al.size(); i++) {
                if(result[j].component.contains(al.getVertex(i))) {
                    *os[j] << ">" << read.id << "\n" << seq << "\n";
                    *alignments[j] << alignment_record;
                    break;
                }
            }
        }
    };
    for(size_t j = 0; j < comps.size(); j++) {
        os[j]->close();
        alignments[j]->close();
        delete os[j];
        delete alignments[j];
        os[j] = nullptr;
        alignments[j] = nullptr;
    }
    return std::move(result);
}

std::vector<Contig> RepeatResolver::ResolveRepeats(logging::Logger &logger, size_t threads,
                              const std::function<bool(const Edge &)> &is_unique) {
    logger.info() << "Splitting dataset and printing subdatasets to disk" << std::endl;
    std::vector<Subdataset> subdatasets = SplitDataset(is_unique);
    logger.info() << "Running repeat resolution" << std::endl;
    std::string COMMAND = "python3 resolution/sequence_graph/path_graph_multik.py -i {} -o {}";
    std::sort(subdatasets.begin(), subdatasets.end());
#pragma omp parallel for default(none) shared(subdatasets, COMMAND)
    for(size_t snum = 0; snum < subdatasets.size(); snum++) {
        Subdataset &subdataset = subdatasets[snum];
        std::experimental::filesystem::path outdir = subdataset.dir / "mltik";
        std::string command = COMMAND;
        command.replace(command.find("{}"), 2, subdataset.dir.string());
        command.replace(command.find("{}"), 2, outdir.string());
        system(command.c_str());
    }
    logger.info() << "Collecting repeat resolution results" << std::endl;
    ParallelRecordCollector<Contig> res(threads);
#pragma omp parallel for default(none) shared(subdatasets, res)
    for(size_t snum = 0; snum < subdatasets.size(); snum++) {
        Subdataset &subdataset = subdatasets[snum];
        std::experimental::filesystem::path outdir = subdataset.dir / "mltik";
        std::experimental::filesystem::path out_fasta;
        bool found = false;
        for (const std::experimental::filesystem::path &f : std::experimental::filesystem::directory_iterator(outdir)) {
            if (endsWith(f.string(), ".fasta")) {
                out_fasta = f;
                found = true;
            }
        }
        io::Library contig_lib = {out_fasta};
        if (!found)
            continue;
        std::vector<Contig> contigs;
        for (StringContig stringContig : io::SeqReader(contig_lib)) {
            Contig contig = stringContig.makeContig();
            if (contig.seq <= !contig.seq) {
                contigs.emplace_back(contig.seq, logging::itos(snum, 1) + "." + contig.id);
            }
        }
        GraphAlignmentStorage storage(dbg);
        for(Contig &contig : contigs) {
            storage.fill(contig);
            res.emplace_back(contig);
        }
        printDot(subdataset.dir / "graph_with_contigs.dot", subdataset.component, storage.labeler());
    }
    logger.info() << "Finished repeat resolution" << std::endl;
    return res.collect();
}

std::vector<Contig> RepeatResolver::CollectResults(logging::Logger &logger, size_t threads, const std::vector<Contig> &contigs,
                                   const std::function<bool(const Edge &)> &is_unique) {
    logger.info() << "Merging results from repeat resolution of subcomponents"<< std::endl;
    logger.info() << "Collecting partial results"<< std::endl;
    ParallelRecordCollector<GraphAlignment> paths(threads);
    omp_set_num_threads(threads);
//#pragma omp parallel default(none) shared(contigs, is_unique, paths)
    for(size_t i = 0; i < contigs.size(); i++) {
        GraphAlignment al = GraphAligner(dbg).align(contigs[i].seq);
        if(al.size() == 1 && is_unique(al[0].contig()))
            continue;
        paths.emplace_back(al);
        paths.emplace_back(al.RC());
    }
    std::vector<GraphAlignment> path_list = paths.collect();
    logger.info() << "Linking contigs"<< std::endl;
    std::unordered_map<dbg::Edge *, size_t> unique_map;
    for(size_t i = 0; i < path_list.size(); i++) {
        Edge & edge = path_list[i].front().contig();
        if(path_list[i].size() == 1 || !is_unique(edge) || path_list[i][0].size() < path_list[i][0].contig().size())
            continue;
        if(unique_map.find(&edge) != unique_map.end()) {
            unique_map[&edge] = size_t(-1);
            unique_map[&edge.rc()] = size_t(-1);
        } else
            unique_map[&edge] = i;
    }
    logger.info() << "Merging contigs"<< std::endl;
    std::vector<Sequence> res;
    std::unordered_set<Edge *> visited_unique;
    for(size_t i = 0; i < path_list.size(); i++) {
        if(!path_list[i].valid()) {
            continue;
        }
        size_t cur = i;
        while(true) {
            Edge &last = path_list[cur].back().contig();
            Segment<Edge> last_seg = path_list[cur].back();
            if(last_seg.size() < last.size() || unique_map.find(&last) == unique_map.end() || unique_map[&last] == size_t(-1))
                break;
            cur = unique_map[&last];
            if(cur == i)
                break;
        }
        cur = cur ^ 1ull;
        size_t start = cur;
        GraphAlignment merged_path;
        while(true) {
            merged_path += path_list[cur].subalignment(1);
            Edge &last = path_list[cur].back().contig();
            Segment<Edge> last_seg = path_list[cur].back();
            path_list[cur] = {};
            if(last_seg.size() < last.size() || unique_map.find(&last) == unique_map.end() || unique_map[&last] == size_t(-1))
                break;
            cur = unique_map[&last];
            if(cur == start)
                break;
        }
        for(Segment<Edge> &seg : merged_path) {
            if(is_unique(seg.contig()) && seg.size() == seg.contig().size()) {
                visited_unique.emplace(&seg.contig());
                visited_unique.emplace(&seg.contig().rc());
            }
        }
        Sequence seq = merged_path.Seq();
        if(!seq < seq)
            seq = !seq;
        res.emplace_back(seq);
        cur = cur ^ 1ull;
        if(!path_list[cur].valid())
            continue;
        start = cur;
        while(true) {
            Edge &last = path_list[cur].back().contig();
            Segment<Edge> last_seg = path_list[cur].back();
            path_list[cur] = {};
            if(last_seg.size() < last.size() || unique_map.find(&last) == unique_map.end() || unique_map[&last] == size_t(-1))
                break;
            cur = unique_map[&last];
            if(cur == start)
                break;
        }
    }
    for(Edge &edge : dbg.edgesUnique()) {
        if(is_unique(edge) && visited_unique.find(&edge) == visited_unique.end()) {
            Sequence seq = edge.start()->seq + edge.seq;
            if(!seq < seq)
                seq = !seq;
            res.emplace_back(seq);
        }
    }
    logger.info() << "Sorting final contigs"<< std::endl;
    std::function<bool(const Sequence &, const Sequence &)> cmp = [](const Sequence &s1, const Sequence &s2){
        if (s1.size() != s2.size())
            return s1.size() > s2.size();
        return s1 < s2;
    };
    std::sort(res.begin(), res.end(), cmp);
    std::vector<Contig> final;
    for(size_t i = 0; i < res.size(); i++) {
        final.emplace_back(res[i], logging::itos(i));
    }
    logger.info() << "Finished clooecting repeat resolution results"<< std::endl;
    return std::move(final);
}

struct RawSeg {
    std::string id;
    size_t left;
    size_t right;

    RawSeg(std::string id, size_t left, size_t right) : id(std::move(id)), left(left), right(right) {}

    bool operator<(const RawSeg &other) const {
        if(id != other.id)
            return id < other.id;
        if(left != other.left)
            return left < other.left;
        return right < other.right;
    }
    bool operator==(const RawSeg &other) const {
        return id == other.id && left == other.left && right == other.right;
    }
};

void PrintFasta(const std::vector<Contig> &contigs, const std::experimental::filesystem::path &path) {
    std::ofstream resos;
    resos.open(path);
    std::vector<Contig> all_contigs;
    for (const Contig &contig : contigs) {
        resos << ">" << contig.id << "\n" << contig.seq << std::endl;
    }
    resos.close();
}

void PrintAlignments(logging::Logger &logger, size_t threads, std::vector<Contig> &contigs,
                     const RecordStorage &readStorage, size_t k, size_t w,
                     const std::experimental::filesystem::path &dir) {
    ensure_dir_existance(dir);
    std::unordered_map<hashing::htype, std::vector<std::pair<Contig *, size_t>>, hashing::alt_hasher<hashing::htype>> position_map;
    std::ofstream refos;
    refos.open(dir/"contigs.fasta");
    std::vector<Contig> all_contigs;
    logger.info() << "Printing compressed assembly to disk" << std::endl;
    for(const Contig & contig : contigs) {
        refos << ">" << contig.id << "\n" << contig.seq << std::endl;
        all_contigs.emplace_back(contig);
        all_contigs.emplace_back(contig.RC());
    }
    refos.close();
    logger.info() << "Aligning reads back to assembly" << std::endl;
    hashing::RollingHash hasher(k, 239);
    for(Contig &contig : all_contigs) {
        for(size_t pos = 1; pos + k <= contig.size(); pos += w) {
            hashing::htype h = hasher.hash(contig.seq, pos);
            position_map[h].emplace_back(&contig, pos);
        }
    }
    ParallelRecordCollector<std::pair<size_t, std::pair<RawSeg, Segment<Contig>>>> result(threads);
    omp_set_num_threads(threads);
#pragma omp parallel for default(none) shared(readStorage, hasher, position_map, k, w, result, std::cout)
    for(size_t i = 0; i < readStorage.size(); i++) {
        const AlignedRead &alignedRead = readStorage[i];
        if(!alignedRead.valid())
            continue;
        Contig read(alignedRead.path.getAlignment().Seq(), alignedRead.id);
        std::vector<std::pair<Segment<Contig>, Segment<Contig>>> res;
        hashing::KWH kwh(hasher, read.seq, 0);
        while (true) {
            if (position_map.find(kwh.fHash()) != position_map.end()) {
                for(std::pair<Contig *, size_t> &pos : position_map[kwh.fHash()]) {
                    size_t shrink_left = 0;
                    if(kwh.pos > pos.second) {
                        shrink_left = kwh.pos - pos.second;
                    }
                    size_t shrink_right = 0;
                    if(-kwh.pos + pos.second + read.size()> pos.first->size())
                        shrink_right = -kwh.pos + pos.second + read.size() - pos.first->size();
                    Segment<Contig> seg_from(read, shrink_left, read.size() - shrink_right);
                    Segment<Contig> seg_to(*pos.first, pos.second - kwh.pos + shrink_left,
                                           pos.second - kwh.pos + read.size() - shrink_right);
                    res.emplace_back(seg_from, seg_to);
                }
            }
            if (!kwh.hasNext())
                break;
            kwh = kwh.next();
        }
        std::sort(res.begin(), res.end());
        res.erase(std::unique(res.begin(), res.end()), res.end());
        for(std::pair<Segment<Contig>, Segment<Contig>> &al : res) {
            if(al.first.seq() == al.second.seq()) {
                RawSeg seg_from(al.first.contig().getId(), al.first.left, al.first.right);
                result.emplace_back(i, std::make_pair(seg_from, al.second));
            }
        }
    }
    std::vector<std::pair<size_t, std::pair<RawSeg, Segment<Contig>>>> final = result.collect();
    __gnu_parallel::sort(final.begin(), final.end());
    logger.info() << "Finished alignment. Printing alignments to " << (dir/"alignments.txt") << std::endl;
    std::ofstream os;
    os.open(dir/"alignments.txt");
    for(auto &rec : final) {
        os << rec.second.first.id << " " << rec.second.second.contig().getId() << " "
           << rec.second.second.left << " " << rec.second.second.right << std::endl;
    }
    os.close();
}
