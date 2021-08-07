#include "subdataset_processing.hpp"

#include <utility>
#include "stdio.h"

std::vector<RepeatResolver::Subdataset> RepeatResolver::SplitDataset(const std::function<bool(const Edge &)> &is_unique) {
    size_t k = dbg.hasher().getK();
    std::vector<Component> comps = ConditionSplitter(is_unique).splitGraph(dbg);
    std::vector<Subdataset> result;
    recreate_dir(dir);
    std::unordered_map<Vertex *, size_t> cmap;
    for(size_t i = 0; i < comps.size(); i++) {
        result.emplace_back(i, comps[i], dir / std::to_string(i));
        for(Vertex &vert : result.back().component.vertices()) {
            cmap[&vert] = result.back().id;
        }
        ensure_dir_existance(result.back().dir);
    }
    for(size_t rnum = 0; rnum < readStorage.size(); rnum++) {
        AlignedRead &read = readStorage[rnum];
        if(!read.valid() || read.path.size() == 1)
            continue;
        GraphAlignment al = read.path.getAlignment();
        result[cmap[&al.getVertex(1)]].reads.emplace_back(rnum);
        for(size_t i = 2; i < al.size(); i++) {
            VERIFY(cmap[&al.getVertex(1)] == cmap[&al.getVertex(i)]);
        }
    }
    return std::move(result);
}

void RepeatResolver::prepareDataset(const RepeatResolver::Subdataset &subdataset) {
    printFasta(subdataset.dir / "graph.fasta", subdataset.component);
    std::ofstream log;
    log.open(subdataset.dir / "dbg.log");
    log << "-k " << dbg.hasher().getK() << std::endl;
    log.close();
    printDot(subdataset.dir / "graph.dot", subdataset.component, readStorage.labeler());
    std::ofstream als;
    als.open(subdataset.dir / "alignments.txt");
    for(size_t rnum : subdataset.reads) {
        AlignedRead &read = readStorage[rnum];
        GraphAlignment al = read.path.getAlignment();
        std::stringstream ss;
        als << read.id << " " << read.path.start().hash() << int(read.path.start().isCanonical())
            << " " << read.path.cpath().str() << "\n";
        CompactPath rc = read.path.RC();
        als  << "-" << read.id << " " << rc.start().hash() << int(rc.start().isCanonical())
            << " " << rc.cpath().str() << "\n";
        std::string alignment_record = ss.str();
    }
    als.close();
}

std::vector<Contig> RepeatResolver::ResolveRepeats(logging::Logger &logger, size_t threads,
                              const std::function<bool(const Edge &)> &is_unique) {
    logger.info() << "Splitting dataset and printing subdatasets to disk" << std::endl;
    std::vector<Subdataset> subdatasets = SplitDataset(is_unique);
    logger.info() << "Running repeat resolution" << std::endl;
    std::string COMMAND = "python3 resolution/sequence_graph/path_graph_multik.py -i {} -o {}";
    std::sort(subdatasets.begin(), subdatasets.end());
#pragma omp parallel for schedule(dynamic, 1) default(none) shared(subdatasets, COMMAND, logger)
    for(size_t snum = 0; snum < subdatasets.size(); snum++) {
        Subdataset &subdataset = subdatasets[snum];
        prepareDataset(subdataset);
        std::experimental::filesystem::path outdir = subdataset.dir / "mltik";
        std::string command = COMMAND;
        command.replace(command.find("{}"), 2, subdataset.dir.string());
        command.replace(command.find("{}"), 2, outdir.string());
        int code = system(command.c_str());
        if(code != 0) {
            logger.info() << "Repeat resolution of component " << subdataset.id << " returned code " << code << std::endl;
            exit(1);
        }
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
                contigs.emplace_back(contig.seq, logging::itos(subdataset.id, 1) + "." + contig.id);
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
                                   const std::experimental::filesystem::path &merging,
                                   const std::function<bool(const Edge &)> &is_unique) {
    logger.info() << "Merging results from repeat resolution of subcomponents"<< std::endl;
    logger.info() << "Collecting partial results"<< std::endl;
    ParallelRecordCollector<AlignedRead> paths(threads);
    omp_set_num_threads(threads);
#pragma omp parallel for default(none) shared(contigs, is_unique, paths)
    for(size_t i = 0; i < contigs.size(); i++) {
        GraphAlignment al = GraphAligner(dbg).align(contigs[i].seq);
        if(al.size() == 1 && is_unique(al[0].contig()))
            continue;
        paths.emplace_back(contigs[i].id, al);
        paths.emplace_back(basic::Reverse(contigs[i].id), al.RC());
    }
    std::vector<AlignedRead> path_list = paths.collect();
    logger.info() << "Linking contigs"<< std::endl;
    std::unordered_map<dbg::Edge *, size_t> unique_map;
    for(size_t i = 0; i < path_list.size(); i++) {
        GraphAlignment al = path_list[i].path.getAlignment();
        Edge & edge = al.front().contig();
        if(path_list[i].path.size() == 1 || !is_unique(edge) || al[0].size() < al[0].contig().size())
            continue;
        if(unique_map.find(&edge) != unique_map.end()) {
            unique_map[&edge] = size_t(-1);
            unique_map[&edge.rc()] = size_t(-1);
        } else {
            unique_map[&edge] = i;
        }
    }
    logger.info() << "Merging contigs"<< std::endl;
    std::vector<Contig> res;
    std::unordered_set<Edge *> visited_unique;
    for(size_t i = 0; i < path_list.size(); i++) {
        if(!path_list[i].valid()) {
            continue;
        }
        size_t cur = i;
        while(true) {
            Segment<Edge> last_seg = path_list[cur].path.getAlignment().back();
            Edge &last = last_seg.contig();
            if(last_seg.size() < last.size() || unique_map.find(&last) == unique_map.end() || unique_map[&last] == size_t(-1))
                break;
            cur = unique_map[&last];
            if(cur == i)
                break;
        }
        cur = cur ^ 1ull;
        size_t start = cur;
        GraphAlignment merged_path = path_list[cur].path.getAlignment().subalignment(0, 1);
        std::vector<std::string> ids;
        size_t clen = merged_path.len();
        ids.emplace_back("(0 -");
        while(true) {
            merged_path += path_list[cur].path.getAlignment().subalignment(1);
            clen += path_list[cur].path.getAlignment().subalignment(1).len();
            ids.emplace_back(path_list[cur].id);
            ids.emplace_back("- " + logging::itos(clen) + ")");
            Segment<Edge> last_seg = path_list[cur].path.getAlignment().back();
            Edge &last = last_seg.contig();
            path_list[cur] = {};
            if(last_seg.size() < last.size() || unique_map.find(&last) == unique_map.end() || unique_map[&last] == size_t(-1))
                break;
            ids.emplace_back("( " + logging::itos(clen - last.size()) + " -");
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
        if(!seq < seq) {
            seq = !seq;
            std::vector<std::string> ids1;
            for(size_t j = 0; j < ids.size(); j++) {
                ids1.emplace_back(basic::Reverse(ids[ids.size() - j - 1]));
            }
            ids = std::move(ids1);
        }
        res.emplace_back(seq, join(" ", ids));
        cur = cur ^ 1ull;
        if(!path_list[cur].valid())
            continue;
        start = cur;
        while(true) {
            Segment<Edge> last_seg = path_list[cur].path.getAlignment().back();
            Edge &last = last_seg.contig();
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
            res.emplace_back(seq, edge.getId());
        }
    }
    logger.info() << "Sorting final contigs"<< std::endl;
    std::function<bool(const Contig &, const Contig &)> cmp = [](const Contig &s1, const Contig &s2){
        if (s1.size() != s2.size())
            return s1.size() > s2.size();
        return s1.seq < s2.seq;
    };
    std::sort(res.begin(), res.end(), cmp);
    std::ofstream merge;
    merge.open(merging);
    std::vector<Contig> final;
    for(size_t i = 0; i < res.size(); i++) {
        final.emplace_back(res[i].seq, logging::itos(i));
        merge << i << "\n" << res[i].id << "\n";
    }
    merge.close();
    logger.info() << "Finished collecting repeat resolution results"<< std::endl;
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
