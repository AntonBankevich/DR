#pragma once
#include "sparse_dbg.hpp"

template<typename htype>
class CompactPath {
private:
    Vertex<htype> *_start;
    Sequence _edges;
    size_t _first_skip;
    size_t _last_skip;
public:
    CompactPath(Vertex<htype> &start, const Sequence& edges, size_t first_skip = 0, size_t last_skip = 0) :
            _start(&start), _edges(edges), _first_skip(first_skip), _last_skip(last_skip) {
    }

    explicit CompactPath(GraphAlignment<htype> &path) :
            _start(&path.getVertex(0)), _first_skip(path.leftSkip()), _last_skip(path.rightSkip()) {
        std::vector<char> edges;
        for(size_t i = 0; i < path.size(); i++) {
            Segment<Edge<htype>> &seg = path[i];
            edges.push_back(seg.contig().seq[0]);
        }
        _edges = Sequence(edges);
    }

    GraphAlignment<htype> getAlignment() {
        std::vector<Segment<Edge<htype>>> path;
        Vertex<htype> *cur = _start;
        for(size_t i = 0; i < _edges.size(); i++) {
            Edge<htype> &edge = cur->getOutgoing(_edges[i]);
            path.emplace_back(edge, 0, edge.size());
            cur = edge.end();
        }
        path.front().left += _first_skip;
        path.back().right += _last_skip;
        return {_start, std::move(path)};
    }

    Path<htype> getPath() {
        std::vector<Edge<htype> *> path;
        Vertex<htype> *cur = _start;
        for(size_t i = 0; i < _edges.size(); i++) {
            Edge<htype> &edge = cur->getOutgoing(_edges[i]);
            path.emplace_back(&edge);
            cur = edge.end();
        }
        return {*_start, std::move(path)};
    }

    Vertex<htype> &start() {
        return *_start;
    }

    const Vertex<htype> &start() const {
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

template<typename htype>
std::ostream& operator<<(std::ostream  &os, const CompactPath<htype> &cpath) {
    return os << cpath.start().hash() << cpath.start().isCanonical() << " " << cpath.cpath() << " " << cpath.leftSkip() << " " << cpath.rightSkip();
}

template<typename htype>
class AlignedRead {
public:
    std::string id;
    CompactPath<htype> path;

    AlignedRead(const Contig &read, GraphAlignment<htype> &_path) : id(read.id), path(_path) {}
    AlignedRead(const Contig &read, CompactPath<htype> &_path) : id(read.id), path(_path) {}
};

template<typename htype>
class AlignedReadStorage {
private:
    SparseDBG<htype> &dbg;
    std::vector<AlignedRead<htype>> reads;
public:
    AlignedReadStorage(SparseDBG<htype> &_dbg) : dbg(_dbg) {
    }

    template<class I>
    void fill(I begin, I end, logging::Logger &logger, size_t threads) {
//        TODO make parallel
        while(begin != end) {
            StringContig read = *begin;
            GraphAlignment<htype> path = dbg.align(read.seq);
            reads.emplace_back(read.makeContig(), path);
            ++begin;
        }
    }
};

template<typename htype>
struct VertexRecord {
private:
    Vertex<htype> &v;
    std::vector<std::pair<Sequence, size_t>> paths;
    size_t zero_cnt = 0;
    size_t cov = 0;
public:
    explicit VertexRecord(Vertex<htype> &_v) : v(_v) {
    }

    VertexRecord(const VertexRecord<htype> &) = delete;

    VertexRecord(VertexRecord<htype> &&other)  noexcept : v(other.v), paths(std::move(other.paths)), zero_cnt(other.zero_cnt), cov(other.cov) {
    }

    VertexRecord<htype> & operator=(const VertexRecord<htype> &) = delete;

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
        for(std::pair<Sequence, size_t> &path : paths) {
            if(path.first == seq) {
                path.second -= 1;
                cov -= 1;
                if(path.second == 0)
                    zero_cnt += 1;
                break;
            }
        }
        if(zero_cnt > paths.size() / 3) {
            std::vector<std::pair<Sequence, size_t>> new_paths;
            for(std::pair<Sequence, size_t> &rec : paths) {
                if(rec.second == 0) {
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

    std::vector<Path<htype>> getAlternatives(const Vertex<htype> &end, double threshold) const {
        lock();
        std::vector<std::pair<Sequence, size_t>> candidates;
        for(const auto & extension : paths) {
            Path<htype> unpacked = CompactPath<htype>(v, extension.first).getPath();
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
        std::vector<Path<htype>> res;
        size_t cnt = 0;
        for(size_t i = 0; i < candidates.size(); i++) {
            if(i > 0 && candidates[i-1].first != candidates[i].first) {
                if(cnt > threshold)
                    res.emplace_back(CompactPath<htype>(v, candidates[i - 1].first).getPath());
                cnt = 0;
            }
            cnt += candidates[i].second;
        }
        if(cnt > threshold)
            res.emplace_back(CompactPath<htype>(v, candidates.back().first).getPath());
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

template<typename htype>
std::ostream& operator<<(std::ostream  &os, const VertexRecord<htype> &rec) {
    return os << rec.str();
}

template<typename htype>
class RecordStorage {
private:
    SparseDBG<htype> &dbg;
    std::vector<AlignedRead<htype>> reads;
    size_t min_len;
    size_t max_len;
    std::unordered_map<const Vertex<htype> *, VertexRecord<htype>> data;

    void processPath(const GraphAlignment<htype> &path, const std::function<void(Vertex<htype> &, const Sequence &)> &task,
                     size_t right = size_t(-1)) {
        CompactPath<htype> cpath(path);
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

    void processPath(CompactPath<htype> &cpath, const std::function<void(Vertex<htype> &, const Sequence &)> &task,
                     size_t right = size_t(-1)) {
        Vertex<htype> *left_vertex = &cpath.start();
        Vertex<htype> *right_vertex = &cpath.start();
        size_t j = 0;
        size_t clen = 0;
        right = std::min(right, cpath.size());
        for (size_t i = 0; i < right; i++) {
            while (j < cpath.size() && clen < max_len) {
                Edge<htype> &edge = right_vertex->getOutgoing(cpath[j]);
                clen += edge.size();
                j++;
                right_vertex = edge.end();
            }
            if (clen < min_len)
                break;
            task(*left_vertex, cpath.cpath().Subseq(i, j));
            Edge<htype> &edge = left_vertex->getOutgoing(cpath[i]);
            clen -= edge.size();
            left_vertex = edge.end();
        }
    }

public:
    RecordStorage(SparseDBG<htype> &_dbg, size_t _min_len, size_t _max_len) : dbg(_dbg), min_len(_min_len), max_len(_max_len) {
        for(auto &it : dbg) {
            data.emplace(&it.second, VertexRecord<htype>(it.second));
            data.emplace(&it.second.rc(), VertexRecord<htype>(it.second.rc()));
        }
    }

    const VertexRecord<htype> &getRecord(const Vertex<htype> &v) const {
        return data.find(&v)->second;
    }

    void addSubpath(const CompactPath<htype> &cpath, size_t right) {
        std::function<void(Vertex<htype> &, const Sequence &)> vertex_task = [this](Vertex<htype> &v, const Sequence &s) {
            data[&v].addPath(s);
        };
        processPath(cpath, vertex_task, right);
    }

    void removeSubpath(const CompactPath<htype> &cpath, size_t right) {
        std::function<void(Vertex<htype> &, const Sequence &)> vertex_task = [this](Vertex<htype> &v, const Sequence &s) {
            data[&v].removePath(s);
        };
        processPath(cpath, vertex_task, right);
    }

    template<class I>
    void fill(I begin, I end, size_t min_read_size, logging::Logger &logger, size_t threads) {
        logger.info() << "Collecting alignments of sequences to the graph" << std::endl;
        std::function<void(Vertex<htype> &, const Sequence &)> vertex_task = [this](Vertex<htype> &v, const Sequence &s) {
            data.find(&v)->second.addPath(s);
        };
        ParallelRecordCollector<AlignedRead<htype>> tmpReads(threads);
        ParallelCounter cnt(threads);
        std::function<void(StringContig &)> read_task = [this, min_read_size, &vertex_task, &tmpReads, &cnt](StringContig & scontig) {
            Contig contig = scontig.makeContig();
            if(contig.size() < min_read_size)
                return;
            GraphAlignment<htype> path = dbg.align(contig.seq);
            GraphAlignment<htype> rcPath = path.RC();
            CompactPath<htype> cpath(path);
            CompactPath<htype> crcPath(rcPath);
            this->processPath(cpath, vertex_task);
            this->processPath(crcPath, vertex_task);
            tmpReads.emplace_back(contig, cpath);
            cnt += cpath.size();
        };
        processRecords(begin, end, logger, threads, read_task);
        reads.insert(reads.end(), tmpReads.begin(), tmpReads.end());
        logger.info() << "Alignment collection finished. Total length of alignments is " << cnt.get() << std::endl;
    }

    typedef typename std::vector<AlignedRead<htype>>::iterator iterator;
    typedef typename std::vector<AlignedRead<htype>>::const_iterator const_iterator;

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
};

