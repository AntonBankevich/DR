#pragma once
#include "sparse_dbg.hpp"
#include <parallel/algorithm>
#include <utility>
using namespace dbg;

class CompactPath {
private:
    Vertex *_start;
    Sequence _edges;
    size_t _first_skip;
    size_t _last_skip;
public:
    CompactPath(Vertex &start, const Sequence& edges, size_t first_skip = 0, size_t last_skip = 0) :
            _start(&start), _edges(edges), _first_skip(first_skip), _last_skip(last_skip) {
    }

    explicit CompactPath(const Path &path, size_t first_skip = 0, size_t last_skip = 0) :
            _start(&path.getVertex(0)), _first_skip(first_skip), _last_skip(last_skip) {
        std::vector<char> edges;
        for(size_t i = 0; i < path.size(); i++) {
            Edge &edge = path[i];
            edges.push_back(edge.seq[0]);
        }
        _edges = Sequence(edges);
    }

    explicit CompactPath(const GraphAlignment &path) :
            _start(&path.getVertex(0)), _first_skip(path.leftSkip()), _last_skip(path.rightSkip()) {
        std::vector<char> edges;
        for(size_t i = 0; i < path.size(); i++) {
            const Segment<Edge> &seg = path[i];
            edges.push_back(seg.contig().seq[0]);
        }
        _edges = Sequence(edges);
    }

    const Sequence &seq() const {
        return _edges;
    }

    GraphAlignment getAlignment() {
        std::vector<Segment<Edge>> path;
        Vertex *cur = _start;
        for(size_t i = 0; i < _edges.size(); i++) {
            Edge &edge = cur->getOutgoing(_edges[i]);
            path.emplace_back(edge, 0, edge.size());
            cur = edge.end();
        }
        path.front().left += _first_skip;
        path.back().right -= _last_skip;
        return {_start, std::move(path)};
    }

    Path getPath() {
        std::vector<Edge *> path;
        Vertex *cur = _start;
        for(size_t i = 0; i < _edges.size(); i++) {
            Edge &edge = cur->getOutgoing(_edges[i]);
            path.emplace_back(&edge);
            cur = edge.end();
        }
        return {*_start, std::move(path)};
    }

    CompactPath RC() {
        return CompactPath(getAlignment().RC());
    }

    Vertex &start() {
        return *_start;
    }

    const Vertex &start() const {
        return *_start;
    }

    const Sequence& cpath() const {
        return _edges;
    }

    size_t leftSkip() const {
        return _first_skip;
    }

    size_t rightSkip() const {
        return _last_skip;
    }

    size_t size() const {
        return _edges.size();
    }

    unsigned char operator[](size_t ind) const {
        return _edges[ind];
    }
};

std::ostream& operator<<(std::ostream  &os, const CompactPath &cpath) {
    return os << cpath.start().hash() << cpath.start().isCanonical() << " " << cpath.cpath() << " " << cpath.leftSkip() << " " << cpath.rightSkip();
}

class AlignedRead {
public:
    std::string id;
    CompactPath path;

    AlignedRead(const Contig &read, GraphAlignment &_path) : id(read.id), path(_path) {}
    AlignedRead(const Contig &read, CompactPath &_path) : id(read.id), path(_path) {}
};

class AlignedReadStorage {
private:
    SparseDBG &dbg;
    std::vector<AlignedRead> reads;
public:
    explicit AlignedReadStorage(SparseDBG &_dbg) : dbg(_dbg) {
    }

    template<class I>
    void fill(I begin, I end, logging::Logger &logger, size_t threads) {
//        TODO make parallel
        while(begin != end) {
            StringContig read = *begin;
            GraphAlignment path = dbg.align(read.makeSequence());
            reads.emplace_back(read.makeContig(), path);
            ++begin;
        }
    }
};

struct VertexRecord {
private:
    typedef std::vector<std::pair<Sequence, size_t>> Storage;
    typedef Storage::const_iterator const_iterator;
    Vertex &v;
    Storage paths;
    size_t zero_cnt = 0;
    size_t cov = 0;
public:
    explicit VertexRecord(Vertex &_v) : v(_v) {
    }

    VertexRecord(const VertexRecord &) = delete;

    VertexRecord(VertexRecord &&other)  noexcept : v(other.v), paths(std::move(other.paths)), zero_cnt(other.zero_cnt), cov(other.cov) {
    }

    const_iterator begin() const {
        return paths.begin();
    }

    const_iterator end() const {
        return paths.end();
    }

    VertexRecord & operator=(const VertexRecord &) = delete;

    void lock() const {
        v.lock();
    }

    void unlock() const {
        v.unlock();
    }

    void addPath(const Sequence &seq) {
        lock();
        cov += 1;
        for(std::pair<Sequence, size_t> &path : paths) {
            if(path.first == seq) {
                path.second += 1;
                if(path.second == 1)
                    zero_cnt -= 1;
                unlock();
                return;
            }
        }
        paths.emplace_back(seq, 1);
        unlock();
    }

    void removePath(const Sequence &seq) {
        lock();
        bool found = false;
        for(std::pair<Sequence, size_t> &path : paths) {
            if(path.first == seq) {
                found = true;
                VERIFY(path.second > 0);
                path.second -= 1;
                cov -= 1;
                if(path.second == 0)
                    zero_cnt += 1;
                break;
            }
        }
        VERIFY(found);
        if(zero_cnt > paths.size() / 3) {
            std::vector<std::pair<Sequence, size_t>> new_paths;
            for(std::pair<Sequence, size_t> &rec : paths) {
                if(rec.second != 0) {
                    new_paths.emplace_back(std::move(rec.first), rec.second);
                }
            }
            std::swap(paths, new_paths);
            zero_cnt = 0;
        }
        unlock();
    }

    size_t countStartsWith(const Sequence &seq) const {
        lock();
        size_t cnt = 0;
        for(const std::pair<Sequence, size_t> &rec : paths) {
            if(rec.first.startsWith(seq)) {
                cnt += rec.second;
            }
        }
        unlock();
        return cnt;
    }

    bool uniqueOut(double threshold) const {
        lock();
        if (paths.empty()) {
            unlock();
            return false;
        }
        size_t best = 0;
        for(size_t i = 1; i < paths.size(); i++) {
            if (paths[i].second > paths[best].second)
                best = i;
        }
        if(paths[best].second < threshold) {
            unlock();
            return false;
        }
        for(const auto & path : paths) {
            if(path.second >= threshold && !paths[best].first.startsWith(path.first) && !path.first.startsWith(paths[best].first)) {
                unlock();
                return false;
            }
        }
        unlock();
        return true;
    }

    std::vector<GraphAlignment> getBulgeAlternatives(const Vertex &end, double threshold) const {
        lock();
        std::vector<std::pair<Sequence, size_t>> candidates;
        for(const auto & extension : paths) {
            Path unpacked = CompactPath(v, extension.first).getPath();
            for(size_t i = 1; i <= unpacked.size(); i++) {
                if(end == unpacked.getVertex(i)){
                    candidates.emplace_back(extension.first.Subseq(0, i), extension.second);
                }
            }
        }
        unlock();
        if(candidates.empty())
            return {};
        std::sort(candidates.begin(), candidates.end());
        std::vector<GraphAlignment> res;
        size_t cnt = 0;
        for(size_t i = 0; i < candidates.size(); i++) {
            if(i > 0 && candidates[i-1].first != candidates[i].first) {
                if(cnt >= threshold)
                    res.emplace_back(CompactPath(v, candidates[i - 1].first).getAlignment());
                cnt = 0;
            }
            cnt += candidates[i].second;
        }
        if(cnt >= threshold)
            res.emplace_back(CompactPath(v, candidates.back().first).getAlignment());
        return std::move(res);
    }

    std::vector<GraphAlignment> getTipAlternatives(size_t len, double threshold) const {
        len += std::max<size_t>(30, len / 20);
        lock();
        std::vector<std::pair<Sequence, size_t>> candidates;
        for(const auto & extension : paths) {
            GraphAlignment unpacked = CompactPath(v, extension.first).getAlignment();
            if(unpacked.len() >= len) {
                unpacked.cutBack(unpacked.len() - len);
                candidates.emplace_back(CompactPath(unpacked).seq(), extension.second);
            }
        }
        unlock();
        if(candidates.empty())
            return {};
        std::sort(candidates.begin(), candidates.end());
        std::vector<GraphAlignment> res;
        size_t cnt = 0;
        for(size_t i = 0; i < candidates.size(); i++) {
            if(i > 0 && candidates[i - 1].first != candidates[i].first) {
                if(cnt >= threshold) {
                    GraphAlignment cp = CompactPath(v, candidates[i - 1].first).getAlignment();
                    cp.cutBack(cp.len() - len);
                    res.emplace_back(cp);
                }
                cnt = 0;
            }
            cnt += candidates[i].second;
        }
        if(cnt >= threshold) {
            GraphAlignment cp = CompactPath(v, candidates.back().first).getAlignment();
            cp.cutBack(cp.len() - len);
            res.emplace_back(cp);
        }
        return std::move(res);
    }

    size_t coverage() const {
        return cov;
    }

    std::string str() const {
        std::stringstream ss;
        lock();
        for(const auto & path : paths) {
            ss << path.first << " " << path.second << std::endl;
        }
        unlock();
        return ss.str();
    }
};

std::ostream& operator<<(std::ostream  &os, const VertexRecord &rec) {
    return os << rec.str();
}

class RecordStorage {
private:
    SparseDBG &dbg;
    std::vector<AlignedRead> reads;
    size_t min_len;
    size_t max_len;
    bool track_cov;
    std::unordered_map<const Vertex *, VertexRecord> data;

    void processPath(const GraphAlignment &path, const std::function<void(Vertex &, const Sequence &)> &task,
                     size_t right = size_t(-1)) {
        CompactPath cpath(path);
        size_t j = 0;
        size_t clen = 0;
        right = std::min(right, cpath.size());
        for (size_t i = 0; i < right; i++) {
            while (j < path.size() && clen < max_len) {
                clen += path[j].contig().size();
                j++;
            }
            if (clen < min_len)
                break;
            task(path.getVertex(i), cpath.cpath().Subseq(i, j));
            clen -= path[i].contig().size();
        }
    }

    void processPath(CompactPath &cpath, const std::function<void(Vertex &, const Sequence &)> &task,
                     const std::function<void(Segment<Edge>)> &edge_task = [](Segment<Edge>){}) {
        Vertex *left_vertex = &cpath.start();
        Vertex *right_vertex = &cpath.start();
        GraphAlignment al = cpath.getAlignment();
        for(size_t i = 0; i < al.size(); i++) {
            Edge &edge = al[i].contig();
            size_t seg_left = i == 0 ? cpath.leftSkip() : 0;
            size_t seg_right = i == cpath.size() - 1 ? edge.size() - cpath.rightSkip() : edge.size();
            edge_task(Segment<Edge>(edge, seg_left, seg_right));
        }
        size_t j = 1;
        size_t clen = 0;
        for (size_t i = 1; i < al.size(); i++) {
            while (j < al.size() && clen < max_len) {
                clen += al[j].contig().size();
                j++;
            }
            if (clen >= min_len)
                task(al.getVertex(i - 1), cpath.cpath().Subseq(i - 1, j));
            clen -= al[i].size();
        }
    }

public:
    RecordStorage(SparseDBG &_dbg, size_t _min_len, size_t _max_len, bool _track_cov = false) :
            dbg(_dbg), min_len(_min_len), max_len(_max_len), track_cov(_track_cov) {
        for(auto &it : dbg) {
            data.emplace(&it.second, VertexRecord(it.second));
            data.emplace(&it.second.rc(), VertexRecord(it.second.rc()));
        }
    }

    const VertexRecord &getRecord(const Vertex &v) const {
        return data.find(&v)->second;
    }

    void addSubpath(CompactPath &cpath) {
        std::function<void(Vertex &, const Sequence &)> vertex_task = [this](Vertex &v, const Sequence &s) {
            data.find(&v)->second.addPath(s);
        };
        std::function<void(Segment<Edge>)> edge_task = [](Segment<Edge> seg){};
        if(track_cov)
            edge_task = [](Segment<Edge> seg){
                seg.contig().incCov(seg.size());
                Edge &edge = seg.contig();
            };
        processPath(cpath, vertex_task, edge_task);
    }

    void removeSubpath(CompactPath &cpath, size_t left = 0, size_t right = size_t(-1)) {
        std::function<void(Vertex &, const Sequence &)> vertex_task = [this](Vertex &v, const Sequence &s) {
            data.find(&v)->second.removePath(s);
        };
        std::function<void(Segment<Edge>)> edge_task = [](Segment<Edge> seg){};
        if(track_cov)
            edge_task = [](Segment<Edge> seg) {
                seg.contig().incCov(size_t(-seg.size()));
            };
        processPath(cpath, vertex_task, edge_task);
    }

    template<class I>
    void fill(I begin, I end, size_t min_read_size, logging::Logger &logger, size_t threads) {
        for(auto & it: dbg) {
            for(Edge &edge : it.second.getOutgoing()) {
                edge.incCov(-edge.intCov());
            }
        }
        logger.info() << "Collecting alignments of sequences to the graph" << std::endl;
        ParallelRecordCollector<AlignedRead> tmpReads(threads);
        ParallelCounter cnt(threads);
        std::function<void(StringContig &)> read_task = [this, min_read_size, &tmpReads, &cnt](StringContig & scontig) {
            Contig contig = scontig.makeContig();
            if(contig.size() < min_read_size)
                return;
            GraphAlignment path = dbg.align(contig.seq);
            GraphAlignment rcPath = path.RC();
            CompactPath cpath(path);
            CompactPath crcPath(rcPath);
            addSubpath(cpath);
            addSubpath(crcPath);
            tmpReads.emplace_back(contig, cpath);
            cnt += cpath.size();
        };
        processRecords(begin, end, logger, threads, read_task);
        reads.insert(reads.end(), tmpReads.begin(), tmpReads.end());
        logger.info() << "Alignment collection finished. Total length of alignments is " << cnt.get() << std::endl;
    }

    void reroute(AlignedRead &alignedRead, const GraphAlignment &initial, const GraphAlignment &corrected) {
        GraphAlignment rcInitial = initial.RC();
        GraphAlignment rcCorrected = corrected.RC();
        CompactPath cInitial(initial);
        CompactPath cCorrected(corrected);
        CompactPath crcInitial(rcInitial);
        CompactPath crcCorrected(rcCorrected);
        this->removeSubpath(cInitial);
        this->removeSubpath(crcInitial);
        this->addSubpath(cCorrected);
        this->addSubpath(crcCorrected);
        alignedRead.path = cCorrected;
    }

    typedef typename std::vector<AlignedRead>::iterator iterator;
    typedef typename std::vector<AlignedRead>::const_iterator const_iterator;

    iterator begin() {
        return reads.begin();
    }

    const_iterator begin() const {
        return reads.begin();
    }

    iterator end() {
        return reads.end();
    }

    const_iterator end() const {
        return reads.end();
    }

    AlignedRead &operator[](size_t ind) {
        return reads[ind];
    }

    const AlignedRead &operator[](size_t ind) const {
        return reads[ind];
    }

    size_t size() const {
        return reads.size();
    }

    void printAlignments(logging::Logger &logger, const std::experimental::filesystem::path &path) const {
        logger.info() << "Printing read to graph alignenments to file " << path << std::endl;
        std::string acgt = "ACGT";
        std::ofstream os;
        os.open(path);
        for(const AlignedRead &read : reads) {
            CompactPath al = read.path;
            os  << read.id << " " << al.start().hash() << int(al.start().isCanonical())
                << " " << al.cpath().str() << "\n";
            al = al.RC();
            os  << "-" << read.id << " " << al.start().hash() << int(al.start().isCanonical())
                << " " << al.cpath().str() << "\n";
        }
        os.close();
    }
};

class CompactPath;

struct GraphError {
    AlignedRead *read;
    size_t from;
    size_t to;
    size_t _size;
    GraphAlignment correction;

    GraphError(AlignedRead *read, size_t from, size_t to, size_t size) :
                    read(read), from(from), to(to), _size(size) {}

    size_t size() const {
        return _size;
    }

    bool operator==(const GraphError &other) const {
        return read == other.read && from == other.from && to == other.to;
    }

};

class ErrorStorage {
private:
    std::vector<GraphError> errors;
    RecordStorage &recs;
public:

    explicit ErrorStorage(RecordStorage &_recs): recs(_recs){
    };

    void fill(logging::Logger &logger, size_t threads, double error_threshold, double extend_threshold) {
        logger.info() << "Correcting low covered regions in reads" << std::endl;
        ParallelRecordCollector<GraphError> result(threads);
#pragma omp parallel for default(none) shared(error_threshold, extend_threshold, result, logger)
        for(size_t read_ind = 0; read_ind < recs.size(); read_ind++) {
            AlignedRead &alignedRead = recs[read_ind];
            std::stringstream ss;
            CompactPath &initial_cpath = alignedRead.path;
            GraphAlignment path = initial_cpath.getAlignment();
            GraphAlignment corrected_path(path.start());
            for(size_t path_pos = 0; path_pos < path.size(); path_pos++) {
                VERIFY_OMP(corrected_path.finish() == path.getVertex(path_pos));
                Edge &edge = path[path_pos].contig();
                if (edge.getCoverage() >= error_threshold) {
                    continue;
                }
                size_t step_back = 0;
                size_t step_front = 0;
                size_t size = edge.size();
                while (step_back < corrected_path.size() &&
                       corrected_path[corrected_path.size() - step_back - 1].contig().getCoverage() <
                       extend_threshold) {
                    size += corrected_path[corrected_path.size() - step_back - 1].size();
                    step_back += 1;
                }
                while (step_front + path_pos + 1 < path.size() &&
                       path[step_front + path_pos + 1].contig().getCoverage() < extend_threshold) {
                    size += path[step_front + path_pos + 1].size();
                    step_front += 1;
                }
                result.emplace_back(&alignedRead, corrected_path.size() - step_back, step_front + path_pos + 1, size);
                path_pos += step_front;
            }
        }
        for(GraphError &error : result)
            errors.emplace_back(error);
    }

    void apply() {
        std::function<bool(const GraphError &, const GraphError &)> compare =
                [](const GraphError &a, const GraphError &b) {
            if (a.read != b.read)
                return a.read < b.read;
            if(a.from != b.from)
                return a.from < b.from;
            return a.to < b.to;
        };
        __gnu_parallel::sort(errors.begin(), errors.end(), compare);
//#pragma omp parallel for default(none) shared(errors, recs)
        for(size_t i = 0; i < errors.size(); i++) {
            if(i > 0 && errors[i].read == errors[i - 1].read) {
                continue;
            }
            std::vector<GraphError> corrected_errors;
            AlignedRead &read = *errors[i].read;
            for(size_t j = i; j < errors.size(); j++) {
                if(errors[j].read != errors[i].read)
                    break;
                if(errors[j].correction.valid()) {
                    corrected_errors.push_back(errors[j]);
                }
            }
            if(corrected_errors.empty())
                continue;
            GraphAlignment initial = read.path.getAlignment();
            std::vector<Segment<Edge>> corrected;
            Vertex &start = corrected_errors[0].from == 0 ? corrected_errors[0].correction.start() : initial.start();
            size_t prev = 0;
            for(GraphError & error : corrected_errors) {
                for(size_t j = prev; j < error.from; j++) {
                    corrected.emplace_back(initial[j]);
                }
                for(size_t j = 0; j < error.size(); j++) {
                    corrected.emplace_back(error.correction[j]);
                }
                prev = error.to;
            }
            for(size_t j = prev; j < initial.size(); j++) {
                corrected.emplace_back(initial[j]);
            }
            GraphAlignment corrected_alignment(&start, std::move(corrected));
            recs.reroute(read, initial, corrected_alignment);
        }
    }

    void sortByLength() {
        std::function<bool(const GraphError &, const GraphError &)> compare =
                [](const GraphError &a, const GraphError &b) {
                    return a._size < b._size;
                };
        __gnu_parallel::sort(errors.begin(), errors.end(), compare);
    }
};

