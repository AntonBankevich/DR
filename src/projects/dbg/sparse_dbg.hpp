//
// Created by anton on 7/22/20.
//

#pragma once
#include "sequences/sequence.hpp"
#include "sequences/seqio.hpp"
#include "common/omp_utils.hpp"
#include "common/logging.hpp"
#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include <common/oneline_utils.hpp>
#include <common/iterator_utils.hpp>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

namespace dbg {
    class Vertex;

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
        std::string getId() const;
        std::string getShortId() const;
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
        void incCov(size_t val) const;
        Sequence kmerSeq(size_t pos) const;
        Sequence suffix(size_t pos) const;
        std::string str() const;
        bool operator==(const Edge &other) const;
        bool operator!=(const Edge &other) const;
        bool operator<(const Edge &other) const;
        bool operator<=(const Edge &other) const;
    };

//    std::ostream& operator<<(std::ostream& os, const Edge& edge);



    class Vertex {
    private:
        friend class SparseDBG;
        mutable std::vector<Edge> outgoing_{};
        Vertex *rc_;
        hashing::htype hash_;
        omp_lock_t writelock;
        size_t coverage_ = 0;
        bool canonical = false;
        bool mark_ = false;
        explicit Vertex(hashing::htype hash, Vertex *_rc);
    public:
        Sequence seq;

        explicit Vertex(hashing::htype hash = 0);
        Vertex(const Vertex &) = delete;
        ~Vertex();

        void mark() {mark_ = true;}
        void unmark() {mark_ = false;}
        bool marked() const {return mark_;}
        hashing::htype hash() const {return hash_;}
        Vertex &rc() {return *rc_;}
        const Vertex &rc() const {return *rc_;}
        void setSequence(const Sequence &_seq);
        void lock() {omp_set_lock(&writelock);}
        void unlock() {omp_unset_lock(&writelock);}
        std::vector<Edge>::iterator begin() const {return outgoing_.begin();}
        std::vector<Edge>::iterator end() const {return outgoing_.end();}
        size_t outDeg() const {return outgoing_.size();}
        size_t inDeg() const {return rc_->outgoing_.size();}
        Edge &operator[](size_t ind) const {return outgoing_[ind];}


        size_t coverage() const;
        bool isCanonical() const;
        bool isCanonical(const Edge &edge) const;
        std::string edgeId(const Edge &edge) const;
        void clear();
        void clearOutgoing();
        void sortOutgoing();
        void checkConsistency() const;
        std::string getId() const;
        std::string getShortId() const;
        void incCoverage();
        void clearSequence();
        Edge &addEdgeLockFree(const Edge &edge);
        void addEdge(const Edge &e);
        Edge &getOutgoing(unsigned char c) const;
        bool hasOutgoing(unsigned char c) const;
        bool isJunction() const;

        bool operator==(const Vertex &other) const;
        bool operator!=(const Vertex &other) const;
        bool operator<(const Vertex &other) const;
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

    struct EdgePosition {
        Edge *edge;
        size_t pos;

        EdgePosition(Edge &_edge, size_t _pos) : edge(&_edge), pos(_pos) {
            VERIFY(pos >= 0 && pos <= edge->size());
        }

        EdgePosition() : edge(nullptr), pos(0) {
        }

        bool isBorder() const {
            return pos == 0 || pos == edge->size();
        }

        std::vector<EdgePosition> step() const {
            if (pos == edge->size()) {
                std::vector<EdgePosition> res;
                Vertex &v = *edge->end();
                for (Edge &next : v) {
                    res.emplace_back(next, 1);
                }
                return std::move(res);
            } else {
                return {{*edge, pos + 1}};
            }
        }

        Sequence kmerSeq() const {
            return edge->kmerSeq(pos);
        }

        unsigned char lastNucl() const {
            return edge->seq[pos - 1];
        }

        EdgePosition RC() const {
            Edge &rc_edge = edge->rc();
            return EdgePosition(rc_edge, edge->size() - pos);
        }
    };

    class SparseDBG {
    public:
        typedef std::unordered_map<hashing::htype , Vertex, hashing::alt_hasher<hashing::htype>> vertex_map_type;
        typedef std::unordered_map<hashing::htype, Vertex, hashing::alt_hasher<hashing::htype>>::iterator vertex_iterator_type;
        typedef std::unordered_map<hashing::htype, EdgePosition, hashing::alt_hasher<hashing::htype>> anchor_map_type;
    private:
//    TODO: replace with perfect hash map? It is parallel, maybe faster and compact.
        vertex_map_type v;
        anchor_map_type anchors;
        hashing::RollingHash hasher_;

//    Be careful since hash does not define vertex. Rc vertices share the same hash
        Vertex &innerAddVertex(hashing::htype h) {
            return v.emplace(std::piecewise_construct, std::forward_as_tuple(h),
                             std::forward_as_tuple(h)).first->second;
        }

    public:

        template<class Iterator>
        SparseDBG(Iterator begin, Iterator end, hashing::RollingHash _hasher) : hasher_(_hasher) {
            while (begin != end) {
                hashing::htype hash = *begin;
                if (v.find(hash) == v.end())
                    addVertex(hash);
                ++begin;
            }
        }

        explicit SparseDBG(hashing::RollingHash _hasher) : hasher_(_hasher) {
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

        SparseDBG SplitGraph(const std::vector<EdgePosition> &breaks) {
            SparseDBG res(hasher_);
            for(auto &it : v) {
                res.addVertex(it.second.seq);
            }
            std::unordered_set<Edge *> broken_edges;
            for(const EdgePosition &epos : breaks) {
                if(!epos.isBorder()) {
                    res.addVertex(epos.kmerSeq());
                    broken_edges.emplace(epos.edge);
                    broken_edges.emplace(&epos.edge->rc());
                }
            }
            for(Edge &edge : edges()) {
                if(broken_edges.find(&edge) == broken_edges.end()) {
                    Vertex &start = res.getVertex(*edge.start());
                    Vertex &end = res.getVertex(*edge.end());
                    Edge new_edge(&start, &end, edge.seq);
                    start.addEdge(new_edge);
                } else {
                    Vertex &newVertex = edge.start()->isCanonical() ? res.getVertex(edge.start()->hash())
                                                                    : res.getVertex(edge.start()->hash()).rc();
                    res.processEdge(newVertex, edge.seq);
                }
            }
            return std::move(res);
        }

        bool containsVertex(const hashing::htype &hash) const {
            return v.find(hash) != v.end();
        }

        void checkConsistency(size_t threads, logging::Logger &logger) {
            logger.info() << "Checking consistency" << std::endl;
            std::function<void(std::pair<const hashing::htype, Vertex> &)> task =
                    [this](std::pair<const hashing::htype, Vertex> &pair) {
                        const Vertex &vert = pair.second;
                        vert.checkConsistency();
                        vert.rc().checkConsistency();
                    };
            processObjects(v.begin(), v.end(), logger, threads, task);
            logger.info() << "Consistency check success" << std::endl;
        }

        void checkSeqFilled(size_t threads, logging::Logger &logger);

        const hashing::RollingHash &hasher() const {
            return hasher_;
        }

        void addVertex(hashing::htype h) {
            innerAddVertex(h);
        }

        Vertex &addVertex(const hashing::KWH &kwh) {
            Vertex &newVertex = innerAddVertex(kwh.hash());
            Vertex &res = kwh.isCanonical() ? newVertex : newVertex.rc();
            res.setSequence(kwh.getSeq());
            return res;
        }

        Vertex &addVertex(const Sequence &seq) {
            return addVertex(hashing::KWH(hasher_, seq, 0));
        }

        Vertex &bindTip(Vertex &start, Edge &tip) {
            Sequence seq = start.seq + tip.seq;
            Vertex &end = addVertex(seq.Subseq(seq.size() - hasher().getK()));
            tip.bindTip(start, end);
            return end;
        }

        Vertex &getVertex(const hashing::KWH &kwh) {
            auto it = v.find(kwh.hash());
            VERIFY(it != v.end());
            if (kwh.isCanonical()) {
                return it->second;
            } else {
                return it->second.rc();
            }
        }

        Vertex &getVertex(const Sequence &seq) {
            return getVertex(hashing::KWH(hasher_, seq, 0));
        }

        Vertex &getVertex(hashing::htype hash) {
            return v.find(hash)->second;
        }

        Vertex &getVertex(const Vertex &other_graph_vertex) {
            if(other_graph_vertex.isCanonical())
                return v.find(other_graph_vertex.hash())->second;
            else
                return v.find(other_graph_vertex.hash())->second.rc();
        }

        std::array<Vertex *, 2> getVertices(hashing::htype hash) {
            Vertex &res = v.find(hash)->second;
            return {&res, &res.rc()};
        }

        const Vertex &getVertex(const hashing::KWH &kwh) const {
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
            ParallelRecordCollector<std::pair<const hashing::htype, EdgePosition>> res(threads);
            std::function<void(Edge &)> task = [&res, w, this](Edge &edge) {
                Vertex &vertex = *edge.start();
                if (edge.size() > w) {
                    Sequence seq = vertex.seq + edge.seq;
//                    Does not run for the first and last kmers.
                    for (hashing::KWH kmer(this->hasher_, seq, 1); kmer.hasNext(); kmer = kmer.next()) {
                        if (kmer.pos % w == 0) {
                            EdgePosition ep(edge, kmer.pos);
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

        void fillAnchors(size_t w, logging::Logger &logger, size_t threads, const std::unordered_set<hashing::htype, hashing::alt_hasher<hashing::htype>> &to_add) {
            logger.info() << "Adding anchors from long edges for alignment" << std::endl;
            ParallelRecordCollector<std::pair<const hashing::htype, EdgePosition>> res(threads);
            std::function<void(Edge &)> task = [&res, w, this, &to_add](Edge &edge) {
                Vertex &vertex = *edge.start();
                if (edge.size() > w || !to_add.empty()) {
                    Sequence seq = vertex.seq + edge.seq;
//                    Does not run for the first and last kmers.
                    for (hashing::KWH kmer(this->hasher_, seq, 1); kmer.hasNext(); kmer = kmer.next()) {
                        if (kmer.pos % w == 0 || to_add.find(kmer.hash()) != to_add.end()) {
                            EdgePosition ep(edge, kmer.pos);
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

        bool isAnchor(hashing::htype hash) const {
            return anchors.find(hash) != anchors.end();
        }

        EdgePosition getAnchor(const hashing::KWH &kwh) {
            if (kwh.isCanonical())
                return anchors.find(kwh.hash())->second;
            else
                return anchors.find(kwh.hash())->second.RC();
        }

        std::vector<hashing::KWH> extractVertexPositions(const Sequence &seq) const {
            std::vector<hashing::KWH> res;
            hashing::KWH kwh(hasher(), seq, 0);
            while (true) {
                if (containsVertex(kwh.hash())) {
                    res.emplace_back(kwh);
                }
                if (!kwh.hasNext())
                    break;
                kwh = kwh.next();
            }
            return std::move(res);
        }

        void printEdge(std::ostream &os, Vertex &start, Edge &edge, bool output_coverage);

        void printDot(std::ostream &os, bool output_coverage) {
            os << "digraph {\nnodesep = 0.5;\n";
            for (std::pair<const hashing::htype, Vertex> &it : this->v) {
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

        void printFastaOld(const std::experimental::filesystem::path &out) {
            std::ofstream os;
            os.open(out);
            for(Edge &edge : edges()) {
                os << ">" << edge.start()->hash() << edge.start()->isCanonical() << "ACGT"[edge.seq[0]] << "\n" <<
                      edge.start()->seq << edge.seq << "\n";
            }
            os.close();
        }

        void processRead(const Sequence &seq) {
            std::vector<hashing::KWH> kmers = extractVertexPositions(seq);
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
            std::vector<hashing::KWH> kmers = extractVertexPositions(seq);
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

        IterableStorage<ApplyingIterator<vertex_iterator_type, Vertex, 2>> vertices(bool unique = false) {
            std::function<std::array<Vertex*, 2>(std::pair<const hashing::htype, Vertex> &)> apply =
                    [unique](std::pair<const hashing::htype, Vertex> &it) -> std::array<Vertex*, 2> {
                if(unique)
                    return {&it.second};
                else
                    return {&it.second, &it.second.rc()};
            };
            ApplyingIterator<vertex_iterator_type, Vertex, 2> begin(v.begin(), v.end(), apply);
            ApplyingIterator<vertex_iterator_type, Vertex, 2> end(v.end(), v.end(), apply);
            return {begin, end};
        }

        IterableStorage<ApplyingIterator<vertex_iterator_type, Vertex, 2>> verticesUnique() {
            return vertices(true);
        }


        IterableStorage<ApplyingIterator<vertex_iterator_type, Edge, 8>> edges(bool unique = false) {
            std::function<std::array<Edge*, 8>(const std::pair<const hashing::htype, Vertex> &)> apply = [unique](const std::pair<const hashing::htype, Vertex> &it) {
                std::array<Edge*, 8> res = {};
                size_t cur = 0;
                for(Edge &edge : it.second) {
                    if(!unique || edge <= edge.rc()) {
                        res[cur] = &edge;
                        cur++;
                    }
                }
                for(Edge &edge : it.second.rc()) {
                    if(!unique || edge <= edge.rc()) {
                        res[cur] = &edge;
                        cur++;
                    }
                }
                return res;
            };
            ApplyingIterator<vertex_iterator_type, Edge, 8> begin(v.begin(), v.end(), apply);
            ApplyingIterator<vertex_iterator_type, Edge, 8> end(v.end(), v.end(), apply);
            return {begin, end};
        }

        IterableStorage<ApplyingIterator<vertex_iterator_type, Edge, 8>> edgesUnique() {
            return edges(true);
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
            std::vector<hashing::htype> todelete;
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
                if (it->second.marked() || (it->second.inDeg() == 0 && it->second.outDeg() == 0)) {
                    it = v.erase(it);
                } else {
                    ++it;
                }
            }
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
