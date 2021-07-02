//
// Created by anton on 7/22/20.
//

#pragma once
#include "sequences/sequence.hpp"
#include "sequences/seqio.hpp"
#include "common/hash_utils.hpp"
#include "common/omp_utils.hpp"
#include "common/logging.hpp"
#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include <vector>
#include <numeric>
#include <unordered_map>
#include <common/oneline_utils.hpp>
#include <unordered_set>

namespace dbg {
    class Vertex;

    class Path;

    class SparseDBG;

    class Edge {
    private:
        Vertex *start_;
        Vertex *end_;
        mutable size_t cov;
    public:
        mutable size_t extraInfo;
        Sequence seq;
        std::string id = "";
        friend class Vertex;
        bool is_reliable = false;
        Edge(Vertex *_start, Vertex *_end, const Sequence &_seq) :
                start_(_start), end_(_end), cov(0), extraInfo(-1), seq(_seq) {
        }
        Vertex *end() const;
        Vertex *start() const;
        size_t getTipSize() const;
        size_t updateTipSize() const;
        void bindTip(Vertex &start, Vertex &end);
        size_t common(const Sequence &other) const;
        size_t size() const;
        double getCoverage() const;
        size_t intCov() const;
        Edge &rc() const;
        Edge &sparseRcEdge() const;
        Path walkForward();
        void incCov(size_t val) const;
        Sequence kmerSeq(size_t pos) const;
        std::string str() const;
        bool operator==(const Edge &other) const;
        bool operator!=(const Edge &other) const;
        bool operator<(const Edge &other) const;
        bool operator<=(const Edge &other) const;
    };

    class Vertex {
    private:
        friend class SparseDBG;

        mutable std::vector<Edge> outgoing_{};
        Vertex *rc_;
        htype hash_;
        omp_lock_t writelock{};
        size_t coverage_ = 0;
        bool canonical = false;
        bool mark_ = false;

        explicit Vertex(htype hash, Vertex *_rc);

    public:
        Sequence seq;

        size_t coverage() const;

        bool isCanonical() const;

        bool isCanonical(const Edge &edge) const;

        std::string edgeId(const Edge &edge) const;

        void mark() {
            mark_ = true;
        }

        void unmark() {
            mark_ = false;
        }

        bool marked() const {
            return mark_;
        }

        void clear() {
            outgoing_.clear();
            rc_->outgoing_.clear();
        }

        void clearOutgoing() {
            outgoing_.clear();
        }

        explicit Vertex(htype hash = 0) : hash_(hash), rc_(new Vertex(hash, this)), canonical(true) {
            omp_init_lock(&writelock);
        }

        void sortOutgoing() {
            std::sort(outgoing_.begin(), outgoing_.end());
        }

        Vertex(const Vertex &) = delete;

        ~Vertex() {
            if (rc_ != nullptr) {
                rc_->rc_ = nullptr;
                delete rc_;
            }
            rc_ = nullptr;
        }

        void checkConsistency() const {
            for (const Edge &edge : outgoing_) {
                if (edge.end() != nullptr) {
                    if (edge.rc().end() != &(this->rc())) {
                        std::cout << this << " " << seq << " " << edge.seq << " " << edge.rc().end() << " "
                                  << &(this->rc()) << std::endl;
                    }
                    VERIFY(edge.rc().end() == &(this->rc()));
                    VERIFY(edge.start_ == this);
                }
            }
        }

        htype hash() const {
            return hash_;
        }

        const Vertex &rc() const {
            return *rc_;
        }

        Vertex &rc() {
            return *rc_;
        }

        void incCoverage() {
#pragma omp atomic update
            coverage_ += 1;
#pragma omp atomic update
            rc().coverage_ += 1;
        }

        Sequence pathSeq(const std::vector<Edge *> &path) const {
            SequenceBuilder sb;
            sb.append(seq);
            for (const Edge *e : path) {
                sb.append(e->seq);
            }
            return sb.BuildSequence();
        }

        void setSequence(const Sequence &_seq) {
            lock();
            if (seq.empty()) {
                if (seq.empty()) {
                    seq = Sequence(_seq.str());
                    unlock();
                    rc_->lock();
                    rc_->seq = !_seq;
                    rc_->unlock();
                } else {
                    unlock();
                }
            } else {
//            if(seq != _seq) {
//                std::cout << seq << std::endl << _seq << std::endl;
//                VERIFY(false);
//            }
//            VERIFY(_seq == seq);
                unlock();
            }
        }

        void clearSequence() {
            if (!seq.empty()) {
                seq = Sequence();
                rc_->seq = Sequence();
            }
        }

        Edge &addEdgeLockFree(const Edge &edge) {
            for (Edge &e : outgoing_) {
                if (edge.size() <= e.size()) {
                    if (edge.seq == e.seq.Subseq(0, edge.size())) {
                        return e;
                    }
                } else if (edge.seq.Subseq(0, e.size()) == e.seq) {
                    e = edge;
                    return e;
                }
            }
            outgoing_.emplace_back(edge);
            return outgoing_.back();
        }

        void addEdge(const Edge &e) {
            omp_set_lock(&writelock);
            addEdgeLockFree(e);
            omp_unset_lock(&writelock);
        }

        void removeEdgesTo(const Vertex &other) {
            lock();
            outgoing_ = oneline::filter<Edge>(outgoing_.begin(), outgoing_.end(),
                                              [&other](Edge &edge) { return edge.end() != &other; });
            unlock();
        }

        void lock() {
            omp_set_lock(&writelock);
        }

        void unlock() {
            omp_unset_lock(&writelock);
        }

        std::vector<Edge>::iterator begin() const {
            return outgoing_.begin();
        }

        std::vector<Edge>::iterator end() const {
            return outgoing_.end();
        }

        Edge &getOutgoing(unsigned char c) const {
            for (Edge &edge : outgoing_) {
                if (edge.seq[0] == c) {
                    return edge;
                }
            }
            std::cout << seq << std::endl;
            std::cout << c << std::endl;
            for (const Edge &edge : outgoing_) {
                std::cout << edge.seq << std::endl;
            }
            VERIFY(false);
            return outgoing_[0];
        }

        bool hasOutgoing(unsigned char c) const {
            for (const Edge &edge : outgoing_) {
                if (edge.seq[0] == c) {
                    return true;
                }
            }
            return false;
        }

        size_t outDeg() const {
            return outgoing_.size();
        }

        size_t inDeg() const {
            return rc_->outgoing_.size();
        }

        bool isJunction() const {
            return outDeg() != 1 || inDeg() != 1;
        }

        Edge &operator[](size_t ind) const {
            return outgoing_[ind];
        }

        bool operator==(const Vertex &other) const {
            return this == &other;
        }

        bool operator!=(const Vertex &other) const {
            return this != &other;
        }

        bool operator<(const Vertex &other) const {
            return hash_ < other.hash_ || (hash_ == other.hash_ && canonical && !other.canonical);
        }

        std::string label() const {
            std::stringstream res;
            if (!isCanonical())
                res << "-";
            res << hash() % 10000000;
            return res.str();
        }
    };

    class EdgePointer {
    private:
        Vertex *vertex;
        unsigned char c;
    public:
        EdgePointer(Edge &edge) : vertex(edge.start()), c(edge.seq[0]) {
        }

        Edge &operator*() const {
            return vertex->getOutgoing(c);
        }

        bool operator==(const EdgePointer &other) {
            return vertex == other.vertex && c == other.c;
        }

        size_t hash() const {
            return std::hash<void*>()(vertex) ^ std::hash<unsigned char>()(c);
        }
    };

    template<class U, class V>
    class PerfectAlignment {
    public:
        Segment<U> seg_from;
        Segment<V> seg_to;

        PerfectAlignment(const Segment<U> &seg_from_, const Segment<V> &seg_to_) : seg_from(seg_from_),
                                                                                   seg_to(seg_to_) {
            VERIFY(seg_from_.size() == seg_to_.size());
        }

        size_t size() {
            return seg_from.size();
        }
    };

    class Path {
    private:
        Vertex *start_;
        std::vector<Edge *> path;
    public:
        Path(Vertex &_start, std::vector<Edge *> _path) : start_(&_start), path(std::move(_path)) {}

        explicit Path(Vertex &_start) : start_(&_start) {}

        Path subPath(size_t from, size_t to) {
            if (from == to)
                return Path(getVertex(from));
            else
                return Path(getVertex(from), std::vector<Edge *>(path.begin() + from, path.begin() + to));
        }

        Path RC() {
            std::vector<Edge *> rcPath;
            for (size_t i = path.size(); i > 0; i--) {
                rcPath.emplace_back(&path[i - 1]->rc());
            }
            return Path(back().end()->rc(), rcPath);
        }

        Edge &operator[](size_t i) const {
            return *path[i];
        }

        Edge &operator[](size_t i) {
            return *path[i];
        }

        Vertex &getVertex(size_t i) {
            VERIFY(i <= path.size());
            if (i == 0)
                return *start_;
            else
                return *path[i - 1]->end();
        }

        Vertex &getVertex(size_t i) const {
            VERIFY(i <= path.size());
            if (i == 0)
                return *start_;
            else
                return *path[i - 1]->end();
        }

        Vertex &start() const {
            return *start_;
        }

        Vertex &finish() const {
            if (path.empty())
                return start();
            else
                return *path.back()->end();
        }

        Edge &back() {
            return *path.back();
        }

        Edge &front() {
            return *path.front();
        }

        double minCoverage() const {
            double res = 100000;
            for (const Edge *edge : path) {
                res = std::min(edge->getCoverage(), res);
            }
            return res;
        }

        Sequence Seq() const {
            return start_->pathSeq(path);
        }

        Sequence truncSeq() const {
            SequenceBuilder sb;
            for (const Edge *e : path) {
                sb.append(e->seq);
            }
            return sb.BuildSequence();
        }

        size_t size() const {
            return path.size();
        }

        typedef typename std::vector<Edge *>::iterator iterator;
        typedef typename std::vector<Edge *>::const_iterator const_iterator;

        iterator begin() {
            return path.begin();
        }

        iterator end() {
            return path.end();
        }

        const_iterator begin() const {
            return path.begin();
        }

        const_iterator end() const {
            return path.end();
        }

        size_t len() const {
            size_t res = 0;
            for (Edge *edge : path)
                res += edge->size();
            return res;
        }

        Path operator+(const Path &other) const {
            VERIFY(finish() == *other.start_);
            std::vector<Edge *> edges = path;
            edges.insert(edges.end(), other.path.begin(), other.path.end());
            return {start(), std::move(edges)};
        }

        void operator+=(Edge &edge) {
            path.emplace_back(&edge);
        }
    };


    class GraphAlignment {
    private:
        Vertex *start_;
        std::vector<Segment<Edge>> als;
    public:
//    TODO change interface
        GraphAlignment(Vertex *_start, std::vector<Segment<Edge>> &&_path) : start_(_start), als(std::move(_path)) {}

        explicit GraphAlignment(Vertex &_start) : start_(&_start) {}

        GraphAlignment() : start_(nullptr) {}

        GraphAlignment RC() const {
            std::vector<Segment<Edge>> path;
            for (size_t i = 0; i < als.size(); i++) {
                Edge &rc_edge = als[i].contig().rc();
                path.emplace_back(rc_edge, rc_edge.size() - als[i].right, rc_edge.size() - als[i].left);
            }
            GraphAlignment res = {&finish().rc(), {path.rbegin(), path.rend()}};
            return res;
        }

        void invalidate() {
            start_ = nullptr;
            als.clear();
        }

        bool valid() const {
            return start_ != nullptr;
        }

        void push_back(const Segment<Edge> &seg) {
            als.push_back(seg);
        }

        void pop_back() {
            als.pop_back();
        }

        void pop_back(size_t len) {
            als.erase(als.end() - len, als.end());
        }

        void cutBack(size_t l) {
            VERIFY(l <= len());
            while (size() > 0 && als.back().size() < l) {
                l -= als.back().size();
                pop_back();
            }
            if (l > 0) {
                als.back().right -= l;
            }
        }


        GraphAlignment subalignment(size_t from, size_t to) {
            return {&getVertex(from), std::vector<Segment<Edge>>(als.begin() + from, als.begin() + to)};
        }

        GraphAlignment &addStep() {
            als.back().right += 1;
            return *this;
        }

        GraphAlignment &addStep(Edge &edge) {
            als.emplace_back(edge, 0, 1);
            return *this;
        }

        GraphAlignment &extend(const Sequence &seq) {
            VERIFY(valid());
            for (size_t cpos = 0; cpos < seq.size(); cpos++) {
                unsigned char c = seq[cpos];
                if (endClosed()) {
                    Vertex &v = finish();
                    if (v.hasOutgoing(c)) {
                        Edge &edge = v.getOutgoing(c);
                        addStep(edge);
                    } else {
                        invalidate();
                        return *this;
                    }
                } else {
                    if (als.back().contig().seq[als.back().right] == c) {
                        addStep();
                    } else {
                        invalidate();
                        return *this;
                    }
                }
            }
            return *this;
        }

        bool endClosed() const {
            return start_ != nullptr && (als.size() == 0 || als.back().right == als.back().contig().size());
        }

        bool startClosed() const {
            return start_ != nullptr && (als.size() == 0 || als.front().left == 0);
        }

        unsigned char lastNucl() const {
            VERIFY(als.back().size() > 0);
            return als.back().contig().seq[als.back().right - 1];
        }

        size_t leftSkip() const {
            return als.size() == 0 ? 0 : als.front().left;
        }

        size_t rightSkip() const {
            return als.size() == 0 ? 0 : als.back().contig().size() - als.back().right;
        }

        std::vector<GraphAlignment> allSteps() {
            if (als.size() != 0 && als.back().right < als.back().contig().size()) {
                GraphAlignment copy = *this;
                return {std::move(copy.addStep())};
            }
            std::vector<GraphAlignment> res;
            Vertex &end = als.size() == 0 ? *start_ : *back().contig().end();
            for (Edge &edge : end) {
                GraphAlignment copy = *this;
                res.emplace_back(std::move(copy.addStep(edge)));
            }
            return res;
        }

        std::vector<GraphAlignment> allExtensions(size_t len) {
            std::vector<GraphAlignment> res = {*this};
            size_t left = 0;
            size_t right = 1;
            for (size_t l = 0; l < len; l++) {
                for (size_t i = left; i < right; i++) {
                    std::vector<GraphAlignment> tmp = res[i].allSteps();
                    res.insert(res.end(), tmp.begin(), tmp.end());
                }
                left = right;
                right = res.size();
            }
            return std::move(res);
        }

        Sequence map(std::unordered_map<const Edge *, Sequence> &edge_map) {
            SequenceBuilder sb;
            bool start = true;
            for (Segment<Edge> &seg : als) {
                auto it = edge_map.find(&seg.contig());
                if (it == edge_map.end()) {
                    if (start) {
                        sb.append((start_->seq + seg.contig().seq).Subseq(seg.left, seg.right + start_->seq.size()));
                        start = false;
                    } else {
                        sb.append(seg.seq());
                    }
                } else {
                    size_t left = start_->seq.size();
                    if (start) {
                        left = 0;
                    }
                    size_t right = start_->seq.size();
                    size_t sz = it->second.size() - start_->seq.size();
                    if (seg.left == 0 && seg.right == seg.contig().size()) {
                        right += sz;
                    } else if (seg.left == 0) {
                        right += std::min(sz, seg.right);
                    } else if (seg.right == seg.contig().size()) {
                        left += sz - std::min(sz, seg.size());
                        right += sz;
                    } else {
                        size_t l = seg.left * sz / seg.contig().size();
                        left += l;
                        right += std::min(l + seg.size(), sz);
                    }
                    sb.append(it->second.Subseq(left, right));
                    start = false;
                }
            }
            return sb.BuildSequence();
        }

        Sequence Seq() const {
            if (als.size() == 0) {
                return {};
            }
            SequenceBuilder sb;
            size_t k = start_->seq.size();
            if (als[0].left >= k)
                sb.append(als[0].contig().seq.Subseq(als[0].left - k, als[0].right));
            else {
                sb.append(start_->seq.Subseq(als[0].left, k));
                sb.append(als[0].contig().seq.Subseq(0, als[0].right));
            }
            for (size_t i = 1; i < als.size(); i++) {
                sb.append(als[i].seq());
            }
            return sb.BuildSequence();
        }

        Sequence truncSeq() const {
            SequenceBuilder sb;
            for (size_t i = 0; i < als.size(); i++) {
                sb.append(als[i].seq());
            }
            return sb.BuildSequence();
        }

        Sequence truncSeq(size_t start_position, size_t size) const {
            SequenceBuilder sb;
            size_t sz = 0;
            for (size_t i = start_position; i < als.size(); i++) {
//            std::cout << i << " " << sz << " " << size << std::endl;
//            std::cout << als[i].contig().size() << " " << als[i].left << " " << als[i].right << " " << als[i].size() << std::endl;
                if (sz + als[i].size() >= size) {
                    sb.append(als[i].seq().Subseq(0, size - sz));
                    break;
                } else {
                    sb.append(als[i].seq());
                    sz += als[i].size();
                }
            }
            return sb.BuildSequence();
        }

        Vertex &start() const {
            return *start_;
        }

        Vertex &finish() const {
            return als.empty() ? *start_ : *als.back().contig().end();
        }

        Segment<Edge> &back() {
            return als.back();
        }

        Segment<Edge> &front() {
            return als.front();
        }

        const Segment<Edge> &operator[](size_t i) const {
            return als[i];
        }

        Segment<Edge> &operator[](size_t i) {
            return als[i];
        }

        Vertex &getVertex(size_t i) const {
            VERIFY(i <= als.size());
            if (i == 0)
                return *start_;
            else
                return *als[i - 1].contig().end();
        }

        Vertex &getVertex(size_t i) {
            VERIFY(i <= als.size());
            if (i == 0)
                return *start_;
            else
                return *als[i - 1].contig().end();
        }

        typename std::vector<Segment<Edge>>::iterator begin() {
            return als.begin();
        }

        typename std::vector<Segment<Edge>>::iterator end() {
            return als.end();
        }

        typename std::vector<Segment<Edge>>::const_iterator begin() const {
            return als.begin();
        }

        typename std::vector<Segment<Edge>>::const_iterator end() const {
            return als.end();
        }

        GraphAlignment subPath(size_t left, size_t right) const {
            if (left == right)
                return GraphAlignment(getVertex(left));
            else
                return {&getVertex(left), {als.begin() + left, als.begin() + right}};
        }

        GraphAlignment reroute(size_t left, size_t right, const GraphAlignment &rerouting) const {
            VERIFY(getVertex(left) == rerouting.start());
            VERIFY(getVertex(right) == rerouting.finish());
            return subPath(0, left) + rerouting + subPath(right, size());
        }

        GraphAlignment reroute(size_t left, size_t right, const Path &rerouting) const {
            VERIFY(getVertex(left) == rerouting.start());
            VERIFY(getVertex(right) == rerouting.finish());
            return subPath(0, left) + rerouting + subPath(right, size());
        }

        void operator+=(const Path &other) {
            if(!valid()) {
                start_ = &other.start();
            }
            VERIFY(finish() == other.getVertex(0));
            for (Edge *edge : other) {
                operator+=(Segment<Edge>(*edge, 0, edge->size()));
            }
        }

        void operator+=(const GraphAlignment &other) {
            if(!valid()) {
                start_ = &other.start();
            }
            VERIFY(finish() == other.getVertex(0));
            for (const Segment<Edge> &al : other) {
                operator+=(al);
            }
        }

        void operator+=(const Segment<Edge> &other) {
            if(!valid()) {
                start_ = other.contig().start();
            }
            if (!als.empty() && als.back().right < als.back().contig().size()) {
                als.back() = als.back() + other;
            } else {
                VERIFY(als.empty() || other.left == 0);
                VERIFY(finish() == *other.contig().start());
                als.push_back(other);
            }
        }

        GraphAlignment operator+(const GraphAlignment &other) const {
            GraphAlignment res = *this;
            res += other;
            return std::move(res);
        }

        GraphAlignment operator+(const Path &other) const {
            GraphAlignment res = *this;
            res += other;
            return std::move(res);
        }

        GraphAlignment operator+(const Segment<Edge> &other) const {
            GraphAlignment res = *this;
            res += other;
            return std::move(res);
        }

        double minCoverage() const {
            double res = 100000;
            for (const Segment<Edge> &seg : als) {
                res = std::min(seg.contig().getCoverage(), res);
            }
            return res;
        }

        Path path() {
            std::vector<Edge *> res;
            for (auto &seg : als) {
                res.push_back(&seg.contig());
            }
            return {*start_, res};
        }

        size_t size() const {
            return als.size();
        }

        size_t len() const {
            size_t res = 0;
            for (auto &seg : als) {
                res += seg.size();
            }
            return res;
        }

        bool operator==(const GraphAlignment &other) const {
            return start_ == other.start_ && als == other.als;
        }

        bool operator!=(const GraphAlignment &other) const {
            return !operator==(other);
        }
    };

    struct EdgePosition {
        Vertex *start;
        Edge *edge;
        size_t pos;

        EdgePosition(Vertex &_start, Edge &_edge, size_t _pos) : start(&_start), edge(&_edge), pos(_pos) {
            VERIFY(pos > 0);
        }

        EdgePosition() : start(nullptr), edge(nullptr), pos(0) {
        }

        GraphAlignment align(const Sequence &seq) {
            GraphAlignment res(start, {{*edge, pos, pos}});
            Edge *cedge = edge;
            size_t epos = pos;
            for (size_t cpos = 0; cpos < seq.size(); cpos++) {
                unsigned char c = seq[cpos];
                if (epos == cedge->size()) {
                    Vertex &v = *cedge->end();
                    if (v.hasOutgoing(c)) {
                        cedge = &v.getOutgoing(c);
                        res.addStep(*cedge);
                        epos = 1;
                    } else {
                        return {};
                    }
                } else {
                    if (cedge->seq[epos] == c) {
                        res.addStep();
                        epos += 1;
                    } else {
                        return {};
                    }
                }
            }
            return std::move(res);
        }

        std::vector<EdgePosition> step() const {
            if (pos == edge->size()) {
                std::vector<EdgePosition> res;
                Vertex &v = *edge->end();
                for (Edge &next : v) {
                    res.emplace_back(v, next, 1);
                }
                return std::move(res);
            } else {
                return {{*start, *edge, pos + 1}};
            }
        }

        unsigned char lastNucl() const {
            return edge->seq[pos - 1];
        }
    };

    class SparseDBG {
    public:
        struct EdgePosition {
            Edge *edge;
            Vertex *start;
            size_t pos;

            EdgePosition(Edge &edge, Vertex &start, size_t pos) : edge(&edge), start(&start), pos(pos) {}

            EdgePosition(const EdgePosition &other) = default;

            EdgePosition RC() const {
                Vertex &s = *start;
                Edge &rc_edge = edge->rc();
                return EdgePosition(rc_edge, edge->end()->rc(), edge->size() - pos);
            }
        };

        typedef std::unordered_map<htype, Vertex, alt_hasher<htype>> vertex_map_type;
        typedef std::unordered_map<htype, EdgePosition, alt_hasher<htype>> anchor_map_type;
    private:
        class EdgeIterator {
            typedef typename vertex_map_type::const_iterator iterator;
            iterator it;
            iterator end;
            bool rc;
            size_t e_num;

            void seek() {
                if(it == end)
                    return;
                {
                    const Vertex &vert = rc ? it->second.rc() : it->second;
                    if (e_num < vert.outDeg()) {
                        return;
                    }
                }
                e_num = 0;
                if(rc) {
                    rc = false;
                    ++it;
                } else {
                    rc = true;
                }
                while(it != end) {
                    const Vertex &vert = rc ? it->second.rc() : it->second;
                    if (e_num < vert.outDeg()) {
                        return;
                    }
                    if(rc) {
                        rc = false;
                        ++it;
                    } else {
                        rc = true;
                    }
                }
            }

        public:
            typedef typename dbg::Edge value_type;

            EdgeIterator(iterator it, iterator end, bool rc, size_t e_num) : it(it), end(end), rc(rc), e_num(e_num) {
                seek();
            }

            Edge &operator*() const {
                Edge *res = nullptr;
                if(rc)
                    res = &it->second.rc()[e_num];
                else
                    res = &it->second[e_num];
                return *res;
            }

            EdgeIterator& operator++() {
                e_num += 1;
                seek();
                return *this;
            }

            EdgeIterator operator++(int) const {
                EdgeIterator other = *this;
                ++other;
                return other;
            }

            bool operator==(const EdgeIterator &other) const {
                return it == other.it && end == other.end && rc == other.rc && e_num == other.e_num;
            }

            bool operator!=(const EdgeIterator &other) const {
                return !operator==(other);
            }
        };

        class EdgeStorage {
        private:
            const SparseDBG &dbg;
        public:
            EdgeStorage(const SparseDBG &dbg) : dbg(dbg) {
            }

            EdgeIterator begin() const {
                return {dbg.begin(), dbg.end(), false, 0};
            }

            EdgeIterator end() const {
                return {dbg.end(), dbg.end(), false, 0};
            }
        };

    private:
//    TODO: replace with perfect hash map? It is parallel, maybe faster and compact.
        vertex_map_type v;
        anchor_map_type anchors;
        RollingHash hasher_;

        std::vector<KWH> extractVertexPositions(const Sequence &seq) const {
            std::vector<KWH> res;
            KWH kwh(hasher_, seq, 0);
            while (true) {
                if (v.find(kwh.hash()) != v.end()) {
                    res.emplace_back(kwh);
                }
                if (!kwh.hasNext())
                    break;
                kwh = kwh.next();
            }
            return std::move(res);
        }

//    Be careful since hash does not define vertex. Rc vertices share the same hash
        Vertex &innerAddVertex(htype h) {
            return v.emplace(std::piecewise_construct, std::forward_as_tuple(h),
                             std::forward_as_tuple(h)).first->second;
        }

    public:

        template<class Iterator>
        SparseDBG(Iterator begin, Iterator end, RollingHash _hasher) : hasher_(_hasher) {
            while (begin != end) {
                htype hash = *begin;
                if (v.find(hash) == v.end())
                    addVertex(hash);
                ++begin;
            }
        }

        explicit SparseDBG(RollingHash _hasher) : hasher_(_hasher) {
        }

        SparseDBG(SparseDBG &&other) = default;

        SparseDBG &operator=(SparseDBG &&other) = default;

//        SparseDBG(SparseDBG &&other) noexcept: hasher_(other.hasher_) {
//            std::swap(v, other.v);
//            std::swap(anchors, other.anchors);
//        }
//
        SparseDBG(const SparseDBG &other) noexcept = delete;

        SparseDBG Subgraph(std::vector<Segment<Edge>> &pieces) {
            SparseDBG res(hasher_);
            for(auto &it : v) {
                res.addVertex(it.second.seq);
            }
            for(Segment<Edge> &seg : pieces) {
                Vertex *left;
                Vertex *right;
                if (seg.left == 0) {
                    left = &res.getVertex(seg.contig().start()->seq);
                } else {
                    left = &res.addVertex(seg.contig().kmerSeq(seg.left));
                }
                Segment<Edge> rcSeg = seg.RC();
                if (rcSeg.left == 0) {
                    right = &res.getVertex(rcSeg.contig().start()->seq);
                } else {
                    right = &res.addVertex(rcSeg.contig().kmerSeq(rcSeg.left));
                }
                left->addEdge(Edge(left, &right->rc(), seg.seq()));
                right->addEdge(Edge(right, &left->rc(), rcSeg.seq()));
            }
            return std::move(res);
        }

        bool containsVertex(const htype &hash) const {
            return v.find(hash) != v.end();
        }

        void checkConsistency(size_t threads, logging::Logger &logger) {
            logger.info() << "Checking consistency" << std::endl;
            std::function<void(std::pair<const htype, Vertex> &)> task =
                    [this](std::pair<const htype, Vertex> &pair) {
                        const Vertex &vert = pair.second;
                        vert.checkConsistency();
                        vert.rc().checkConsistency();
                    };
            processObjects(v.begin(), v.end(), logger, threads, task);
            logger.info() << "Consistency check success" << std::endl;
        }

        void checkSeqFilled(size_t threads, logging::Logger &logger) {
            logger.info() << "Checking vertex sequences" << std::endl;
            std::function<void(std::pair<const htype, Vertex> &)> task =
                    [&logger](std::pair<const htype, Vertex> &pair) {
                        const Vertex &vert = pair.second;
                        if (vert.seq.empty() || vert.rc().seq.empty()) {
                            logger.info() << "Sequence not filled " << pair.first << std::endl;
                            VERIFY(false);
                        }
                        if (!vert.isCanonical()) {
                            logger.info() << "Canonical vertex marked not canonical " << pair.first << std::endl;
                            VERIFY(false);
                        }
                        if (vert.rc().isCanonical()) {
                            logger.info() << "Noncanonical vertex marked canonical " << pair.first << std::endl;
                            VERIFY(false);
                        }
                    };
            processObjects(v.begin(), v.end(), logger, threads, task);
            logger.info() << "Vertex sequence check success" << std::endl;
        }

        const RollingHash &hasher() const {
            return hasher_;
        }

        void addVertex(htype h) {
            innerAddVertex(h);
        }

        Vertex &addVertex(const KWH &kwh) {
            Vertex &newVertex = innerAddVertex(kwh.hash());
            Vertex &res = kwh.isCanonical() ? newVertex : newVertex.rc();
            res.setSequence(kwh.getSeq());
            return res;
        }

        Vertex &addVertex(const Sequence &seq) {
            return addVertex(KWH(hasher_, seq, 0));
        }

        Vertex &bindTip(Vertex &start, Edge &tip) {
            Sequence seq = start.seq + tip.seq;
            Vertex &end = addVertex(seq.Subseq(seq.size() - hasher().getK()));
            tip.bindTip(start, end);
            return end;
        }

        Vertex &getVertex(const KWH &kwh) {
            auto it = v.find(kwh.hash());
            VERIFY(it != v.end());
            if (kwh.isCanonical()) {
                return it->second;
            } else {
                return it->second.rc();
            }
        }

        Vertex &getVertex(const Sequence &seq) {
            return getVertex(KWH(hasher_, seq, 0));
        }

        Vertex &getVertex(htype hash) {
            return v.find(hash)->second;
        }

        std::vector<Vertex *> getVertices(htype hash) {
            Vertex &res = v.find(hash)->second;
            return {&res, &res.rc()};
        }

        const Vertex &getVertex(const KWH &kwh) const {
            auto it = v.find(kwh.hash());
            VERIFY(it != v.end());
            if (kwh.isCanonical()) {
                return it->second;
            } else {
                return it->second.rc();
            }
        }

        void fillAnchors(size_t w, logging::Logger &logger, size_t threads) {
            logger.info() << "Adding anchors from long edges for alignment" << std::endl;
            ParallelRecordCollector<std::pair<const htype, EdgePosition>> res(threads);
            std::function<void(Edge &)> task = [&res, w, this](Edge &edge) {
                Vertex &vertex = *edge.start();
                if (edge.size() > w) {
                    Sequence seq = vertex.seq + edge.seq;
//                    Does not run for the first and last kmers.
                    for (KWH kmer(this->hasher_, seq, 1); kmer.hasNext(); kmer = kmer.next()) {
                        if (kmer.pos % w == 0) {
                            EdgePosition ep(edge, vertex, kmer.pos);
                            if (kmer.isCanonical())
                                res.emplace_back(kmer.hash(), ep);
                            else {
                                res.emplace_back(kmer.hash(), ep.RC());
                            }
                        }
                    }
                }
            };
            EdgeStorage edgeIt = edges();
            processObjects(edgeIt.begin(), edgeIt.end(), logger, threads, task);
            for (auto &tmp : res) {
                anchors.emplace(tmp);
            }
            logger.info() << "Added " << anchors.size() << " anchors" << std::endl;
        }

        void fillAnchors(size_t w, logging::Logger &logger, size_t threads, const std::unordered_set<htype, alt_hasher<htype>> &to_add) {
            logger.info() << "Adding anchors from long edges for alignment" << std::endl;
            ParallelRecordCollector<std::pair<const htype, EdgePosition>> res(threads);
            std::function<void(Edge &)> task = [&res, w, this, to_add](Edge &edge) {
                Vertex &vertex = *edge.start();
                if (edge.size() > w || !to_add.empty()) {
                    Sequence seq = vertex.seq + edge.seq;
//                    Does not run for the first and last kmers.
                    for (KWH kmer(this->hasher_, seq, 1); kmer.hasNext(); kmer = kmer.next()) {
                        if (kmer.pos % w == 0 || to_add.find(kmer.hash()) != to_add.end()) {
                            EdgePosition ep(edge, vertex, kmer.pos);
                            if (kmer.isCanonical())
                                res.emplace_back(kmer.hash(), ep);
                            else {
                                res.emplace_back(kmer.hash(), ep.RC());
                            }
                        }
                    }
                }
            };
            processObjects(edges().begin(), edges().end(), logger, threads, task);
            for (auto &tmp : res) {
                anchors.emplace(tmp);
            }
            logger.info() << "Added " << anchors.size() << " anchors" << std::endl;
        }

        bool isAnchor(htype hash) const {
            return anchors.find(hash) != anchors.end();
        }

        EdgePosition getAnchor(const KWH &kwh) {
            if (kwh.isCanonical())
                return anchors.find(kwh.hash())->second;
            else
                return anchors.find(kwh.hash())->second.RC();
        }

        GraphAlignment align(const Sequence &seq, Edge *edge_to, size_t pos_to) {
            size_t k = hasher().getK();
            size_t cur = k;
            GraphAlignment res;
            while(cur < seq.size()) {
                size_t len = std::min(seq.size() - cur, edge_to->size() - pos_to);
                res += Segment<Edge>(*edge_to, pos_to, pos_to + len);
                cur += len;
                if(cur < seq.size()) {
                    edge_to = &edge_to->end()->getOutgoing(seq[cur]);
                    pos_to = 0;
                }
            }
            return res;
        }

        GraphAlignment align(const Sequence &seq) {
            std::vector<KWH> kmers = extractVertexPositions(seq);
            std::vector<Segment<Edge>> res;
            if (kmers.size() == 0) {
                KWH kwh(hasher_, seq, 0);
                while (true) {
                    if (isAnchor(kwh.hash())) {
                        EdgePosition pos = getAnchor(kwh);
                        VERIFY(kwh.pos < pos.pos);
                        VERIFY(pos.pos + seq.size() - kwh.pos <= pos.edge->size() + hasher_.getK());
                        Segment<Edge> seg(*pos.edge, pos.pos - kwh.pos, pos.pos + seq.size() - kwh.pos - hasher_.getK());
                        return {pos.start, std::vector<Segment<Edge>>({seg})};
                    }
                    if (!kwh.hasNext()) {
#pragma omp critical
                        {
                            std::cout << "Error: could not align sequence " << seq.size() << std::endl;
                            std::cout << seq << std::endl;
                            abort();
                        };
                        return {nullptr, std::move(res)};
                    }
                    kwh = kwh.next();
                }
            }
            Vertex *prestart = &getVertex(kmers.front());
            if (kmers.front().pos > 0) {
                Vertex &rcstart = prestart->rc();
                if (!rcstart.hasOutgoing(seq[kmers.front().pos - 1] ^ 3)) {
                    std::cout << "No outgoing for start" << std::endl << seq << std::endl <<
                              kmers.front().pos << " " << seq[kmers.front().pos - 1] << std::endl
                              << kmers.front().getSeq() << std::endl;
                    VERIFY(false);
                }
                Edge &rcedge = rcstart.getOutgoing(seq[kmers.front().pos - 1] ^ 3);
                prestart = &rcedge.end()->rc();
                Edge &edge = rcedge.rc();
                VERIFY(edge.size() >= kmers.front().pos);
                Segment<Edge> seg(edge, edge.size() - kmers.front().pos, edge.size());
                res.emplace_back(seg);
            }
            for (const KWH &kmer : kmers) {
                if (kmer.pos + hasher_.getK() < seq.size()) {
                    Vertex &vertex = getVertex(kmer);
                    if (!vertex.hasOutgoing(seq[kmer.pos + hasher_.getK()])) {
                        std::cout << "No outgoing for middle" << std::endl << seq << std::endl <<
                                  kmer.pos << " " << size_t(seq[kmer.pos + hasher_.getK()]) << std::endl
                                  << kmer.getSeq() << std::endl;
                        std::cout << vertex.hash() << " " << vertex.outDeg() << " " << vertex.inDeg() << std::endl;
                        for (const Edge &e : vertex) {
                            std::cout << e.seq << std::endl;
                        }
                        VERIFY(false);
                    }
                    Edge &edge = vertex.getOutgoing(seq[kmer.pos + hasher_.getK()]);
                    Segment<Edge> seg(edge, 0, std::min(seq.size() - kmer.pos - hasher_.getK(), edge.size()));
                    res.emplace_back(seg);
                }
            }
            return {prestart, std::move(res)};
        }

        std::vector<PerfectAlignment<Contig, Edge>> carefulAlign(Contig &contig) {
            Sequence seq = contig.seq;
            std::vector<PerfectAlignment<Contig, Edge>> res;
            KWH kwh(hasher_, seq, 0);
            size_t k = hasher_.getK();
            while (true) {
                if (res.empty() || kwh.pos >= res.back().seg_from.right) {
                    if (containsVertex(kwh.hash())) {
                        Vertex &vertex = getVertex(kwh);
                        Vertex &rcVertex = vertex.rc();
                        if ((res.empty() || kwh.pos > res.back().seg_from.right)
                            && kwh.pos > 0 && rcVertex.hasOutgoing(seq[kwh.pos - 1] ^ 3)) {
                            Edge &edge = rcVertex.getOutgoing(seq[kwh.pos - 1] ^ 3);
                            size_t len = 1;
                            while (len < edge.size() && len < kwh.pos && edge.seq[len] == (seq[kwh.pos - len - 1] ^ 3))
                                len += 1;
                            res.emplace_back(Segment<Contig>(contig, kwh.pos - len, kwh.pos),
                                             Segment<Edge>(edge.rc(), edge.size() - len, edge.size()));
                        }
                        if (kwh.pos + k < seq.size() && vertex.hasOutgoing(seq[kwh.pos + k])) {
                            Edge &edge = vertex.getOutgoing(seq[kwh.pos + k]);
                            size_t len = 1;
                            while (len < edge.size() && kwh.pos + k + len < seq.size() &&
                                   edge.seq[len] == seq[kwh.pos + k + len])
                                len += 1;
                            res.emplace_back(Segment<Contig>(contig, kwh.pos, kwh.pos + len),
                                             Segment<Edge>(edge, 0, len));
                        }
                    } else if ((res.empty() || kwh.pos > res.back().seg_from.right) && isAnchor(kwh.hash())) {
                        typename SparseDBG::EdgePosition pos = getAnchor(kwh);
//                TODO replace this code with a call to expand method of PerfectAlignment class after each edge is marked by its full sequence
                        Edge &edge = *pos.edge;
                        Vertex &start = *pos.start;
                        CompositeSequence edge_seq({start.seq, edge.seq});
                        size_t left_from = kwh.pos;
                        size_t right_from = kwh.pos + k;
                        size_t left_to = pos.pos;
                        size_t right_to = pos.pos + k;
                        while (left_from > 0 && left_to > 0 && edge_seq[left_to - 1] == seq[left_from - 1]) {
                            left_from -= 1;
                            left_to -= 1;
                        }
                        while (right_from < seq.size() && right_to < edge_seq.size() &&
                               seq[right_from] == edge_seq[right_to]) {
                            right_from += 1;
                            right_to += 1;
                        }
                        if (left_to - left_from > k) {
                            res.emplace_back(Segment<Contig>(contig, left_from, right_from - k),
                                             Segment<Edge>(edge, left_to, right_to - k));
                        }
                    }
                }
                if (!kwh.hasNext())
                    break;
                kwh = kwh.next();
            }
            return std::move(res);
        }

        std::vector<PerfectAlignment<Edge, Edge>> oldEdgeAlign(Edge &contig) {
            Sequence seq = contig.start()->seq + contig.seq;
            std::vector<PerfectAlignment<Edge, Edge>> res;
            KWH kwh(hasher_, seq, 0);
            size_t k = hasher_.getK();
            while (true) {
                if(!kwh.hasNext())
                    break;
                if (res.empty() || kwh.pos >= res.back().seg_from.right) {
                    Edge *edge = nullptr;
                    size_t pos = 0;
                    if (containsVertex(kwh.hash())) {
                        Vertex &start = getVertex(kwh);
                        if(start.hasOutgoing(seq[kwh.pos + k]))
                            edge = &getVertex(kwh).getOutgoing(seq[kwh.pos + k]);
                    }
                    if (edge == nullptr && isAnchor(kwh.hash())) {
                        typename SparseDBG::EdgePosition gpos = getAnchor(kwh);
                        edge = gpos.edge;
                        pos = gpos.pos;
                    }
                    if(edge != nullptr) {
                        size_t len = std::min(contig.size() - kwh.pos, edge->size() - pos);
                        res.emplace_back(Segment<Edge>(contig, kwh.pos, kwh.pos + len), Segment<Edge>(*edge, pos, pos + len));
                    }
                }
                kwh = kwh.next();
            }
            return std::move(res);
        }

        void printEdge(std::ostream &os, Vertex &start, Edge &edge, bool output_coverage) {
            Vertex &end = *edge.end();
            os << "\"";
            if (!start.isCanonical())
                os << "-";
            os << start.hash() % 100000 << "\" -> \"";
            if (!end.isCanonical())
                os << "-";
            if (output_coverage)
                os << end.hash() % 100000 << "\" [label=\"" << edge.size() << "(" << edge.getCoverage() << ")\"]\n";
            else
                os << end.hash() % 100000 << "\" [label=\"" << edge.size() << "\"]\n";
        }

        void printDot(std::ostream &os, bool output_coverage) {
            os << "digraph {\nnodesep = 0.5;\n";
            for (std::pair<const htype, Vertex> &it : this->v) {
                Vertex &start = it.second;
                for (Edge &edge : start) {
                    Vertex &end = *edge.end();
                    printEdge(os, start, edge, output_coverage);
                    if (v.find(end.hash()) == v.end()) {
                        printEdge(os, end.rc(), edge.rc(), output_coverage);
                    }
                }
                for (Edge &edge : start.rc()) {
                    Vertex &end = *edge.end();
                    printEdge(os, start.rc(), edge, output_coverage);
                    if (v.find(end.hash()) == v.end()) {
                        printEdge(os, end.rc(), edge.rc(), output_coverage);
                    }
                }
            }
            os << "}\n";
        }


        void processRead(const Sequence &seq) {
            std::vector<KWH> kmers = extractVertexPositions(seq);
            if (kmers.size() == 0) {
                std::cout << seq << std::endl;
            }
            VERIFY(kmers.size() > 0);
            std::vector<Vertex *> vertices;
            for (size_t i = 0; i < kmers.size(); i++) {
                vertices.emplace_back(&getVertex(kmers[i]));
                if (i == 0 || vertices[i] != vertices[i - 1]) {
                    vertices.back()->setSequence(kmers[i].getSeq());
                    vertices.back()->incCoverage();
                }
            }
            for (size_t i = 0; i + 1 < vertices.size(); i++) {
//            TODO: if too memory heavy save only some of the labels
                VERIFY(kmers[i].pos + hasher_.getK() <= seq.size())
                if (i > 0 && vertices[i] == vertices[i - 1] && vertices[i] == vertices[i + 1] &&
                    (kmers[i].pos - kmers[i - 1].pos == kmers[i + 1].pos - kmers[i].pos) &&
                    kmers[i + 1].pos - kmers[i].pos < hasher_.getK()) {
                    continue;
                }
                vertices[i]->addEdge(Edge(vertices[i], vertices[i + 1], Sequence(seq.Subseq(kmers[i].pos + hasher_.getK(),
                                                                                            kmers[i + 1].pos +
                                                                                            hasher_.getK()).str())));
                vertices[i + 1]->rc().addEdge(Edge(&vertices[i + 1]->rc(), &vertices[i]->rc(),
                                                   !Sequence(seq.Subseq(kmers[i].pos, kmers[i + 1].pos).str())));
            }
            if (kmers.front().pos > 0) {
                vertices.front()->rc().addEdge(Edge(&vertices.front()->rc(), nullptr, !(seq.Subseq(0, kmers[0].pos))));
            }
            if (kmers.back().pos + hasher_.getK() < seq.size()) {
                vertices.back()->addEdge(
                        Edge(vertices.back(), nullptr, seq.Subseq(kmers.back().pos + hasher_.getK(), seq.size())));
            }
        }

        void processEdge(Vertex &vertex, Sequence old_seq) {
            Sequence seq = vertex.seq + old_seq;
            std::vector<KWH> kmers = extractVertexPositions(seq);
            VERIFY(kmers.front().pos == 0 && kmers.back().pos == old_seq.size());
            std::vector<Vertex *> vertices;
            for (size_t i = 0; i < kmers.size(); i++) {
                vertices.emplace_back(&getVertex(kmers[i]));
            }
            for (size_t i = 0; i + 1 < vertices.size(); i++) {
//            TODO: if too memory heavy save only some of the labels
                VERIFY(kmers[i].pos + hasher_.getK() <= seq.size())
                if (i > 0 && vertices[i] == vertices[i - 1] && vertices[i] == vertices[i + 1] &&
                    (kmers[i].pos - kmers[i - 1].pos == kmers[i + 1].pos - kmers[i].pos) &&
                    kmers[i + 1].pos - kmers[i].pos < hasher_.getK()) {
                    continue;
                }
                vertices[i]->addEdge(
                        Edge(vertices[i], vertices[i + 1], old_seq.Subseq(kmers[i].pos, kmers[i + 1].pos)));
            }
        }

        size_t size() const {
            return v.size();
        }

        EdgeStorage edges() const {
            return {*this};
        }

        void printCoverageStats(logging::Logger &logger) const {
            std::vector<size_t> cov(100);
            std::vector<size_t> cov_tips(100);
            std::vector<std::vector<size_t>> cov_ldist(100);
            std::vector<std::vector<size_t>> cov_ldist_tips(100);
            for (size_t i = 0; i < cov_ldist.size(); i++) {
                cov_ldist[i].resize(30);
            }
            for (size_t i = 0; i < cov_ldist.size(); i++) {
                cov_ldist_tips[i].resize(30);
            }
            std::vector<size_t> covLen(100);
            std::vector<size_t> covLen_tips(100);
            for (auto &val: v) {
                const Vertex &tmp = val.second;
                for (const Edge &edge : tmp) {
                    cov[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                    covLen[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
                    cov_ldist[std::min(size_t(edge.getCoverage()), cov.size() - 1)][std::min(edge.size() / 50,
                                                                                             cov_ldist[0].size() -
                                                                                             1)] += 1;
                    if (edge.getTipSize() < hasher_.getK() * 2) {
                        cov_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                        covLen_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
                        cov_ldist_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)][std::min(edge.size() / 50,
                                                                                                      cov_ldist[0].size() -
                                                                                                      1)] += 1;
                    }
                }
                for (const Edge &edge : tmp.rc()) {
                    cov[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                    covLen[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
                    cov_ldist[std::min(size_t(edge.getCoverage()), cov.size() - 1)][std::min(edge.size() / 50,
                                                                                             cov_ldist[0].size() -
                                                                                             1)] += 1;
                    if (edge.getTipSize() < hasher_.getK() * 2) {
                        cov_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                        covLen_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
                        cov_ldist_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)][std::min(edge.size() / 50,
                                                                                                      cov_ldist[0].size() -
                                                                                                      1)] += 1;
                    }
                }
            }
//        logger.info() << "Distribution of coverages:" << std::endl;
//        for(size_t i = 0; i < cov.size(); i++) {
//            logger << i << " " << cov[i] << " " << covLen[i];
//            for(size_t val : cov_ldist[i]) {
//                logger << " " << val;
//            }
//            logger << std::endl;
//            logger << i << " " << cov_tips[i] << " " << covLen_tips[i];
//            for(size_t val : cov_ldist_tips[i]) {
//                logger << " " << val;
//            }
//            logger << std::endl;
//        }
        }

        typename vertex_map_type::iterator begin() {
            return v.begin();
        }

        typename vertex_map_type::iterator end() {
            return v.end();
        }

        typename vertex_map_type::const_iterator begin() const {
            return v.begin();
        }

        typename vertex_map_type::const_iterator end() const {
            return v.end();
        }

        void removeIsolated() {
            vertex_map_type newv;
            std::vector<htype> todelete;
            for (auto it = v.begin(); it != v.end();) {
                if (it->second.outDeg() == 0 && it->second.inDeg() == 0) {
                    it = v.erase(it);
                } else {
                    ++it;
                }
            }
        }

        void removeMarked() {
            for (auto it = v.begin(); it != v.end();) {
                if (it->second.marked()) {
                    it = v.erase(it);
                } else {
                    ++it;
                }
            }
        }

        void printFasta(std::ostream &out) const {
            size_t cnt = 0;
            for (const auto &it : v) {
                const Vertex &vertex = it.second;
                VERIFY(!vertex.seq.empty());
                for (size_t i = 0; i < vertex.outDeg(); i++) {
                    const Edge &edge = *(vertex.begin() + i);
                    Sequence tmp = vertex.seq + edge.seq;
                    Vertex &end = *edge.end();
                    out << ">" << cnt << "_" << vertex.hash() << int(vertex.isCanonical()) <<
                        "_" << end.hash() << int(end.isCanonical()) << "_" << edge.size() << "_" << edge.getCoverage()
                        << std::endl;
                    cnt++;
                    out << tmp.str() << "\n";
                }
                const Vertex &rcvertex = vertex.rc();
                for (size_t i = 0; i < rcvertex.outDeg(); i++) {
                    const Edge &edge = rcvertex[i];
                    Sequence tmp = rcvertex.seq + edge.seq;
                    Vertex &end = *edge.end();
                    out << ">" << cnt << "_" << rcvertex.hash() << int(rcvertex.isCanonical()) <<
                        "_" << end.hash() << int(end.isCanonical()) << "_" << edge.size() << "_" << edge.getCoverage()
                        << std::endl;
                    cnt++;
                    out << tmp.str() << "\n";
                }
            }
        }

        void printGFA(std::ostream &out, bool calculate_coverage) const {
            out << "H\tVN:Z:1.0" << std::endl;
            size_t cnt = 0;
            for (const auto &it : v) {
                for (const Vertex *pv : {&it.second, &it.second.rc()}) {
                    const Vertex &vertex = *pv;
                    VERIFY(!vertex.seq.empty());
                    for (const Edge &edge : vertex) {
                        if (vertex.isCanonical(edge)) {
                            if (calculate_coverage)
                                out << "S\t" << vertex.edgeId(edge) << "\t" << vertex.seq << edge.seq
                                    << "\tKC:i:" << edge.intCov() << std::endl;
                            else
                                out << "S\t" << vertex.edgeId(edge) << "\t" << vertex.seq << edge.seq << std::endl;
                        }
                    }
                }
            }
            for (const auto &it : v) {
                const Vertex &vertex = it.second;
                for (const Edge &out_edge : vertex) {
                    std::string outid = vertex.edgeId(out_edge);
                    bool outsign = vertex.isCanonical(out_edge);
                    for (const Edge &inc_edge : vertex.rc()) {
                        std::string incid = vertex.rc().edgeId(inc_edge);
                        bool incsign = !vertex.rc().isCanonical(inc_edge);
                        out << "L\t" << incid << "\t" << (incsign ? "+" : "-") << "\t" << outid << "\t"
                            << (outsign ? "+" : "-") << "\t" << hasher_.getK() << "M" << std::endl;
                    }
                }
            }
        }

        template<class Iterator>
        void fillSparseDBGEdges(Iterator begin, Iterator end, logging::Logger &logger, size_t threads,
                                const size_t min_read_size) {
            typedef typename Iterator::value_type ContigType;
            logger.info() << "Starting to fill edges" << std::endl;
            std::function<void(ContigType &)> task = [this, min_read_size](ContigType &contig) {
                Sequence seq = contig.makeSequence();
                if (seq.size() >= min_read_size)
                    processRead(seq);
            };
            processRecords(begin, end, logger, threads, task);
            logger.info() << "Sparse graph edges filled." << std::endl;
        }

        template<class Iterator>
        void refillSparseDBGEdges(Iterator begin, Iterator end, logging::Logger &logger, size_t threads) {
            logger.info() << "Starting to fill edges" << std::endl;
            std::function<void(std::pair<Vertex *, Sequence> &)> task = [this](std::pair<Vertex *, Sequence> &contig) {
                processEdge(*contig.first, contig.second);
            };
            processObjects(begin, end, logger, threads, task);
            logger.info() << "Sparse graph edges filled." << std::endl;
        }

        static SparseDBG
        loadDBGFromFasta(const io::Library &lib, RollingHash &hasher, logging::Logger &logger, size_t threads) {
            logger.info() << "Loading graph from fasta" << std::endl;
            io::SeqReader reader(lib);
            ParallelRecordCollector<Sequence> sequences(threads);
            ParallelRecordCollector<htype> vertices(threads);
            std::function<void(StringContig &)> collect_task = [&sequences, &vertices, hasher](StringContig &contig) {
                Sequence seq = contig.makeSequence();
                KWH start(hasher, seq, 0);
                KWH end(hasher, !seq, 0);
                vertices.add(start.hash());
                vertices.add(end.hash());
                sequences.add(seq);
            };
            processRecords(reader.begin(), reader.end(), logger, threads, collect_task);
            SparseDBG res(vertices.begin(), vertices.end(), hasher);
            reader.reset();
            res.fillSparseDBGEdges(sequences.begin(), sequences.end(), logger, threads, hasher.getK() + 1);
            logger.info() << "Finished loading graph" << std::endl;
            return std::move(res);
        }
    };
}

namespace std {
    template <>
    struct hash<dbg::EdgePointer> {
        std::size_t operator()(const dbg::EdgePointer &p) const {
            return p.hash();
        }
    };
}
