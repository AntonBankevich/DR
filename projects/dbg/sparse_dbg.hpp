//
// Created by anton on 7/22/20.
//

#pragma once
#include "sequences/sequence.hpp"
#include "sequences/seqio.hpp"
#include "common/omp_utils.hpp"
#include "common/logging.hpp"
#include "rolling_hash.hpp"
#include "hash_utils.hpp"
#include <vector>
#include <numeric>
#include <unordered_map>
#include <common/oneline_utils.hpp>
#include <unordered_set>

template<class htype>
class Vertex;

template<class htype>
class Edge {
private:
    Vertex<htype> *end_;
    mutable size_t cov;
public:
    mutable size_t extraInfo;
    Sequence seq;
    std::string id;
    friend class Vertex<htype>;

    Edge(Vertex<htype> *_end, const Sequence &_seq) :
            end_(_end), cov(0), extraInfo(-1), seq(_seq), id("") {
    }

    Vertex<htype> *end() const {
        return end_;
    }

    size_t getTipSize() const {
        return extraInfo;
    }

    size_t updateTipSize() const {
        size_t new_val = 0;
        if(extraInfo == size_t(-1) && end_->inDeg() == 1) {
            for (const Edge<htype> & other : end_->getOutgoing()) {
                other.end_->lock();
                new_val = std::max(new_val, other.extraInfo);
                other.end_->unlock();
            }
            if(new_val != size_t(-1))
                new_val += size();
            end_->lock();
            extraInfo = new_val;
            end_->unlock();
        }
        return new_val;
    }

    void bindTip(Vertex<htype> & start, Vertex<htype> & end) {
        VERIFY(end_ == nullptr);
        end_ = &end;
        Sequence rcseq = !(start.seq + seq);
        end.rc().addEdgeLockFree(Edge<htype>(&start.rc(), rcseq.Subseq(start.seq.size())));
    }

    size_t common(Sequence other) const {
        size_t res = 0;
        while(res < seq.size() && res < other.size() && seq[res] == other[res]) {
            res += 1;
        }
        return res;
    }

    size_t size() const {
        return seq.size();
    }

    double getCoverage() const {
        return double(cov) / size();
    }

    double intCov() const {
        return cov;
    }


    void incCov(size_t val) const {
#pragma omp atomic
        cov += val;
    }

    bool operator==(const Edge<htype> &other) const {
        return this == &other;
    }

    bool operator<(const Edge<htype> &other) const {
        return this->seq < other.seq;
    }
};

template<typename htype>
class SparseDBG;

template<class htype>
class Vertex {
private:
    friend class SparseDBG<htype>;
    std::vector<Edge<htype>> outgoing_{};
    Vertex * rc_;
    htype hash_;
    omp_lock_t writelock{};
    size_t coverage_ = 0;
    bool canonical = false;
    bool mark_ = false;


    explicit Vertex(htype hash, Vertex *_rc) : hash_(hash), rc_(_rc), canonical(false) {
        omp_init_lock(&writelock);
    }

public:
    Sequence seq;

    size_t converage() const {
        return coverage_;
    }

    bool isCanonical() const {
        return canonical;
    }

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

    explicit Vertex(htype hash = 0) : hash_(hash), rc_(new Vertex<htype>(hash, this)), canonical(true) {
        omp_init_lock(&writelock);
    }

    Vertex(Vertex &&other) noexcept : rc_(other.rc_), hash_(other.hash_), canonical(other.canonical) {
        std::swap(outgoing_, other.outgoing_);
        std::swap(seq, other.seq);
        omp_init_lock(&writelock);
        if(other.rc_ != nullptr) {
            other.rc_->rc_ = this;
            other.rc_ = nullptr;
        }
        for(Edge<htype> &edge : rc_->outgoing_) {
            if(edge.end() == &other) {
                edge.end_ = this;
            } else if (edge.end() != nullptr){
                for(Edge<htype> &back : edge.end()->rc().outgoing_) {
                    if(back.end() == &other) {
                        back.end_ = this;
                    }
                }
            }
        }
    }

    void sortOutgoing() {
        std::sort(outgoing_.begin(), outgoing_.end());
    }

    Vertex(const Vertex &) = delete;

    ~Vertex() {
        if(rc_ != nullptr) {
            rc_->rc_ = nullptr;
            delete rc_;
        }
        rc_ = nullptr;
    }

    void checkConsistency() const {
        for(const Edge<htype> & edge : outgoing_) {
            if(edge.end() != nullptr) {
                if(this->rcEdge(edge).end() != &(this->rc())) {
                    std::cout << this << " " << seq << " " << edge.seq << " " << rcEdge(edge).end() << " " << &(this->rc()) <<std::endl;
                }
                VERIFY(this->rcEdge(edge).end() == &(this->rc()));
            }
        }
    }

    htype hash() const {
        return hash_;
    }

    const Vertex<htype> & rc() const {
        return *rc_;
    }

    Vertex<htype> & rc() {
        return *rc_;
    }

    void incCoverage() {
#pragma omp atomic update
        coverage_ += 1;
#pragma omp atomic update
        rc().coverage_ += 1;
    }

    Edge<htype>& rcEdge(const Edge<htype> & edge) {
        Vertex &vend = edge.end()->rc();
        char c;
        if(edge.size() > seq.size()) {
            c = (!edge.seq)[seq.size()];
        } else {
            c = (!seq)[seq.size() - edge.size()];
        }
        return vend.getOutgoing(c);
    }

    Edge<htype>& sparseRcEdge(const Edge<htype> & edge) {
        Vertex &vend = edge.end()->rc();
        VERIFY(seq.size() > 0);
        for(Edge<htype> &candidate : vend.getOutgoing()) {
            if(candidate.end() == rc_ && candidate.size() == edge.size() &&
                    (edge.size() <= seq.size() || candidate.seq.startsWith((!edge.seq).Subseq(seq.size())))) {
                return candidate;
            }
        }
        std::cout << seq + edge.seq << std::endl;
        std::cout << seq << std::endl;
        std::cout << vend.seq << std::endl;
        for(Edge<htype> &candidate : vend.getOutgoing()) {
            if(candidate.end() == rc_ && candidate.size() == edge.size() &&
               (edge.size() <= seq.size() || candidate.seq.startsWith((!edge.seq).Subseq(seq.size())))) {
                std::cout << vend.seq + candidate.seq << std::endl;
            }
        }
        VERIFY(false);
        return vend.getOutgoing()[0];
    }

    const Edge<htype>& sparseRcEdge(const Edge<htype> & edge) const {
        Vertex &vend = edge.end()->rc();
        VERIFY(seq.size() > 0);
        for(Edge<htype> &candidate : vend.getOutgoing()) {
            if(candidate.end() == rc_ && candidate.size() == edge.size() &&
               (edge.size() <= seq.size() || candidate.seq.startsWith((!edge.seq).Subseq(seq.size())))) {
                return candidate;
            }
        }
        std::cout << seq + edge.seq << std::endl;
        std::cout << seq << std::endl;
        std::cout << vend.seq << std::endl;
        for(Edge<htype> &candidate : vend.getOutgoing()) {
            if(candidate.end() == rc_ && candidate.size() == edge.size() &&
               (edge.size() <= seq.size() || candidate.seq.startsWith((!edge.seq).Subseq(seq.size())))) {
                std::cout << vend.seq + candidate.seq << std::endl;
            }
        }
        VERIFY(false);
        return vend.getOutgoing()[0];
    }

    const Edge<htype>& rcEdge(const Edge<htype> & edge) const {
        const Vertex &vend = edge.end()->rc();
        char c;
        if(edge.size() > seq.size()) {
            c = (!edge.seq)[seq.size()];
        } else {
            c = (!seq)[seq.size() - edge.size()];
        }
        return vend.getOutgoing(c);
    }

    Sequence pathSeq(const std::vector<Edge<htype>> & path) const {
        SequenceBuilder sb;
        sb.append(seq);
        for(const Edge<htype> &e : path) {
            sb.append(e.seq);
        }
        return sb.BuildSequence();
    }

    std::vector<Edge<htype>> walkForward(const Edge<htype> & edge) const {
        std::vector<Edge<htype>> res;
        res.push_back(edge);
        Vertex<htype> *next = edge.end();
        while(next != nullptr && next != this && !next->isJunction()) {
            res.push_back(next->getOutgoing()[0]);
            next = next->getOutgoing()[0].end();
        }
        return res;
    }

    void setSequence(const Sequence &_seq) {
        lock();
        if(seq.empty()) {
            if(seq.empty()) {
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

    void addEdgeLockFree(const Edge<htype> &edge) {
        for(Edge<htype> & e : outgoing_) {
            if (edge.size() <= e.size()) {
                if (edge.seq == e.seq.Subseq(0, edge.size())) {
                    return;
                }
            } else if (edge.seq.Subseq(0, e.size()) == e.seq) {
                e = edge;
                return;
            }
        }
        outgoing_.emplace_back(edge);
    }

    void addEdge(const Edge<htype> &e) {
        omp_set_lock(&writelock);
        addEdgeLockFree(e);
        omp_unset_lock(&writelock);
    }

    void removeEdgesTo(const Vertex<htype> & other) {
        lock();
        outgoing_ = oneline::filter(other.begin(), other.end(), [&](const Edge<htype> & edge) {return edge.end() != &other;});
        unlock();
    }

    void lock() {
        omp_set_lock(&writelock);
    }

    void unlock() {
        omp_unset_lock(&writelock);
    }

    const std::vector<Edge<htype>> &getOutgoing() const {
        return outgoing_;
    }

    const Edge<htype> & getOutgoing(char c) const {
        for(const Edge<htype> &edge : outgoing_) {
            if(edge.seq[0] == c) {
                return edge;
            }
        }
        std::cout << seq << std::endl;
        std::cout << c << std::endl;
        for(const Edge<htype> &edge : outgoing_) {
            std::cout << edge.seq << std::endl;
        }
        VERIFY(false);
        return getOutgoing()[0];
    }

    Edge<htype> & getOutgoing(char c) {
        for(Edge<htype> &edge : outgoing_) {
            if(edge.seq[0] == c) {
                return edge;
            }
        }
        std::cout << seq << std::endl;
        std::cout << size_t(c) << std::endl;
        for(Edge<htype> &edge : outgoing_) {
            std::cout << edge.seq << std::endl;
        }
        VERIFY(false);
        return getOutgoing()[0];
    }

    bool hasOutgoing(char c) const {
        for(const Edge<htype> &edge : outgoing_) {
            if(edge.seq[0] == c) {
                return true;
            }
        }
        return false;
    }

    std::vector<Edge<htype>> &getOutgoing() {
        return outgoing_;
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

    bool operator==(const Vertex<htype> & other) const {
        return this == &other;
    }

    bool operator!=(const Vertex<htype> & other) const {
        return this != &other;
    }

    bool operator<(const Vertex<htype> & other) const {
        return hash_ < other.hash_ || (hash_ == other.hash_ && canonical && !other.canonical);
    }
};


template<typename htype>
class Path {
private:
    Vertex<htype> *start_;
    std::vector<Edge<htype>> path;
public:
    Path(Vertex<htype> &_start, std::vector<Edge<htype>> _path) : start_(&_start), path(std::move(_path)){}

    Path<htype> subPath(size_t from, size_t to) {
        return Path(getVertex(from), std::vector<Edge<htype>>(path.begin() + from, path.begin() + to));
    }

    Path<htype> RC() {
        std::vector<Edge<htype>> rcPath;
        for(size_t i = path.size(); i > 0; i--) {
            rcPath.emplace_back(getVertex(i - 1).rcEdge(path[i - 1]));
        }
        return Path(back().end()->rc(), rcPath);
    }

    const Edge<htype> &operator[](size_t i) const {
        return path[i];
    }

    Edge<htype> &operator[](size_t i) {
        return path[i];
    }

    Vertex<htype> &getVertex(size_t i) {
        VERIFY(i <= path.size());
        if(i == 0)
            return *start_;
        else
            return *path[i - 1].end();
    }

    const Vertex<htype> &getVertex(size_t i) const {
        VERIFY(i <= path.size());
        if(i == 0)
            return *start_;
        else
            return *path[i - 1].end();
    }

    Vertex<htype> &start() {
        return *start_;
    }

    Edge<htype> &back() {
        return path.back();
    }

    Sequence Seq() const {
        return start_->pathSeq(path);
    }

    size_t size() const {
        return path.size();
    }
};

template<typename htype>
class GraphAlignment {
private:
    Vertex<htype> *start;
    std::vector<Segment<Edge<htype>>> als;
public:
    GraphAlignment(Vertex<htype> *_start, std::vector<Segment<Edge<htype>>> _path) : start(_start), als(std::move(_path)){}

    const Segment<Edge<htype>> &operator[](size_t i) const {
        return als[i];
    }

    typename std::vector<Segment<Edge<htype>>>::iterator begin() {
        return als.begin();
    }

    typename std::vector<Segment<Edge<htype>>>::iterator end() {
        return als.end();
    }

    typename std::vector<Segment<Edge<htype>>>::const_iterator begin() const {
        return als.begin();
    }

    typename std::vector<Segment<Edge<htype>>>::const_iterator end() const {
        return als.end();
    }

    Path<htype> path() {
        std::vector<Edge<htype>> res;
        for(auto & seg : als) {
            res.push_back(seg.contig());
        }
        return {*start, res};
    }

    size_t size() const {
        return als.size();
    }
};

template<typename htype>
class SparseDBG {
public:
    struct EdgePosition {
        Edge<htype> *edge;
        Vertex<htype> *start;
        size_t pos;

        EdgePosition(Edge<htype> &edge, Vertex<htype> &start, size_t pos) : edge(&edge), start(&start), pos(pos) {}
        EdgePosition(const EdgePosition &other) : edge(other.edge), start(other.start), pos(other.pos) {}

        EdgePosition RC() {
            Vertex<htype> &s = *start;
            Edge<htype> &rc_edge = s.rcEdge(*edge);
            return EdgePosition(rc_edge, edge->end()->rc(), edge->size() - pos);
        }
    };
private:
//    TODO: replace with perfect hash map? It is parallel, maybe faster and compact.
    std::unordered_map<htype, Vertex<htype>> v;
    std::unordered_map<htype, EdgePosition> anchors;
    const RollingHash<htype> hasher_;

    std::vector<KWH<htype>> extractVertexPositions(const Sequence &seq) const {
        std::vector<KWH<htype>> res;
        KWH<htype> kwh(hasher_, seq, 0);
        while(true) {
            if(v.find(kwh.hash()) != v.end()) {
                res.emplace_back(kwh);
            }
            if (!kwh.hasNext())
                break;
            kwh = kwh.next();
        }
        return std::move(res);
    }

//    Be careful since hash does not define vertex. Rc vertices share the same hash
    Vertex<htype> & innerAddVertex(htype h) {
        return v.emplace(h, Vertex<htype>(h)).first->second;
    }

public:

    template<class Iterator>
    SparseDBG(Iterator begin, Iterator end, RollingHash<htype> _hasher) : hasher_(_hasher) {
        while(begin != end) {
            htype hash = *begin;
            if(v.find(hash) == v.end())
                addVertex(hash);
            ++begin;
        }
    }

    SparseDBG(RollingHash<htype> _hasher) : hasher_(_hasher) {
    }

    SparseDBG(SparseDBG<htype> &&other) noexcept = default;

    SparseDBG &operator=(SparseDBG<htype> &&other)   noexcept {
        std::swap(v, other.v);
        std::swap(anchors, other.anchors);
        hasher = other.hasher;
    }

    bool containsVertex(const htype &hash) const {
        return v.find(hash) != v.end();
    }

    void checkConsistency(size_t threads, logging::Logger &logger) {
        std::cout << "Checking consistency" << std::endl;
        std::function<void(std::pair<const htype, Vertex<htype>> &)> task =
                [this](std::pair<const htype, Vertex<htype>> & pair) {
                    const Vertex<htype> & vert = pair.second;
                    vert.checkConsistency();
                    vert.rc().checkConsistency();
                };
        processObjects(v.begin(), v.end(), logger, threads, task);
        std::cout << "Consistency check success" << std::endl;
    }

    void checkSeqFilled(size_t threads, logging::Logger &logger) {
        std::cout << "Checking vertex sequences" << std::endl;
        std::function<void(std::pair<const htype, Vertex<htype>> &)> task =
                [](std::pair<const htype, Vertex<htype>> & pair) {
                    const Vertex<htype> & vert = pair.second;
                    if(vert.seq.empty() || vert.rc().seq.empty()) {
                        std::cout << "Sequence not filled " << pair.first << std::endl;
                        VERIFY(false);
                    }
                    if(!vert.isCanonical()) {
                        std::cout << "Canonical vertex marked not canonical " << pair.first << std::endl;
                        VERIFY(false);
                    }
                    if(vert.rc().isCanonical()) {
                        std::cout << "Noncanonical vertex marked canonical " << pair.first << std::endl;
                        VERIFY(false);
                    }
                };
        processObjects(v.begin(), v.end(), logger, threads, task);
        std::cout << "Vertex sequence check success" << std::endl;
    }

    const RollingHash<htype> &hasher() const {
        return hasher_;
    }

    void addVertex(htype h) {
        innerAddVertex(h);
    }

    Vertex<htype> &addVertex(const KWH<htype> &kwh) {
        Vertex<htype> &newVertex = innerAddVertex(kwh.hash());
        Vertex<htype> &res = kwh.isCanonical() ? newVertex : newVertex.rc();
        res.setSequence(kwh.getSeq());
        return res;
    }

    Vertex<htype> &addVertex(const Sequence &seq) {
        return addVertex(KWH<htype>(hasher_, seq, 0));
    }

    Vertex<htype> &bindTip(Vertex<htype> & start, Edge<htype> &tip) {
        Sequence seq = start.seq + tip.seq;
        Vertex<htype> &end = addVertex(seq.Subseq(seq.size() - hasher().k));
        tip.bindTip(start, end);
        return end;
    }

    Vertex<htype> &getVertex(const KWH<htype> &kwh) {
        VERIFY(v.find(kwh.hash()) != v.end());
        if(kwh.isCanonical()) {
            return v[kwh.hash()];
        } else {
            return v[kwh.hash()].rc();
        }
    }

    Vertex<htype> &getVertex(htype hash) {
        return v.find(hash)->second;
    }

    const Vertex<htype> &getVertex(const KWH<htype> &kwh) const {
        VERIFY(v.find(kwh.hash()) != v.end());
        if(kwh.isCanonical()) {
            return v.find(kwh.hash())->second;
        } else {
            return v.find(kwh.hash())->second.rc();
        }
    }

    void fillAnchors(size_t w, logging::Logger &logger, size_t threads) {
        logger << "Adding anchors from long edges for alignment" << std::endl;
        ParallelRecordCollector<std::pair<const htype, EdgePosition>> res(threads);
        std::function<void(std::pair<const htype, Vertex<htype>> &)> task = [&res, w, this] (std::pair<const htype, Vertex<htype>> &iter) {
            Vertex<htype> &vertex = iter.second;
            for(Edge<htype> &edge : vertex.getOutgoing()) {
                if(edge.size() > w) {
                    Sequence seq = vertex.seq + edge.seq;
//                    Does not run for the first and last kmers.
                    for(KWH<htype> kmer(this->hasher_, seq, 1); kmer.hasNext(); kmer = kmer.next()) {
                        if (kmer.pos % w == 0) {
                            EdgePosition ep(edge, vertex, kmer.pos);
                            if(kmer.isCanonical())
                                res.emplace_back(kmer.hash(), ep);
                            else {
                                res.emplace_back(kmer.hash(), ep.RC());
                            }
                        }
                    }
                }
            }
            for(Edge<htype> &edge : vertex.rc().getOutgoing()) {
                if(edge.size() > w) {
                    Sequence seq = vertex.rc().seq + edge.seq;
//                    Does not run for the first and last kmers.
                    for(KWH<htype> kmer(this->hasher_, seq, 1); kmer.hasNext(); kmer = kmer.next()) {
                        if (kmer.pos % w == 0) {
                            EdgePosition ep(edge, vertex.rc(), kmer.pos);
                            if(kmer.isCanonical())
                                res.emplace_back(kmer.hash(), ep);
                            else {
                                res.emplace_back(kmer.hash(), ep.RC());
                            }
                        }
                    }
                }
            }
        };
        processObjects(begin(), end(), logger, threads, task);
        for(auto & tmp : res) {
            anchors.emplace(tmp);
        }
        logger << "Added " << anchors.size() << " anchors" << std::endl;
    }

    bool isAnchor(htype hash) const {
        return anchors.find(hash) != anchors.end();
    }

    EdgePosition getAnchor(const KWH<htype> &kwh) {
        if(kwh.isCanonical())
            return anchors.find(kwh.hash())->second;
        else
            return anchors.find(kwh.hash())->second.RC();
    }

    GraphAlignment<htype> align(const Sequence & seq) {
        std::vector<KWH<htype>> kmers = extractVertexPositions(seq);
        std::vector<Segment<Edge<htype>>> res;
        if(kmers.size() == 0) {
            KWH<htype> kwh(hasher_, seq, 0);
            while(true) {
                if(isAnchor(kwh.hash())) {
                    EdgePosition pos = getAnchor(kwh);
                    VERIFY(kwh.pos < pos.pos);
                    VERIFY(pos.pos + seq.size() - kwh.pos <= pos.edge->size() + hasher_.k);
                    Segment<Edge<htype>> seg(*pos.edge, pos.pos - kwh.pos, pos.pos + seq.size() - kwh.pos - hasher_.k);
                    return {pos.start, std::vector<Segment<Edge<htype>>>({seg})};
                }
                if (!kwh.hasNext()) {
#pragma omp critical
                    {
                        std::cout << "Gopa" << std::endl;
                        std::cout << seq << std::endl;
                    };
                    return {nullptr, res};
                }
                kwh = kwh.next();
            }
        }
        Vertex<htype> *prestart = &getVertex(kmers.front());
        if (kmers.front().pos > 0) {
            const Vertex<htype> &rcstart = prestart->rc();
            if(!rcstart.hasOutgoing(seq[kmers.front().pos - 1] ^ 3)) {
                std::cout << "No outgoing for start" << std::endl << seq << std::endl <<
                        kmers.front().pos << " " << seq[kmers.front().pos - 1] << std::endl
                        << kmers.front().getSeq() << std::endl;
                VERIFY(false);
            }
            const Edge<htype> &rcedge = rcstart.getOutgoing(seq[kmers.front().pos - 1] ^ 3);
            prestart = &rcedge.end()->rc();
            const Edge<htype> &edge = rcstart.rcEdge(rcedge);
            VERIFY(edge.size() >= kmers.front().pos);
            Segment<Edge<htype>> seg(edge, edge.size() - kmers.front().pos, edge.size());
            res.emplace_back(seg);
        }
        for(const KWH<htype> & kmer : kmers) {
            if (kmer.pos + hasher_.k < seq.size()) {
                const Vertex<htype> &vertex = getVertex(kmer);
                if(!vertex.hasOutgoing(seq[kmer.pos + hasher_.k])) {
                    std::cout << "No outgoing for middle" << std::endl << seq << std::endl <<
                              kmer.pos << " " << size_t(seq[kmer.pos + hasher_.k]) << std::endl
                              << kmer.getSeq() << std::endl;
                    std::cout << vertex.hash() << " " << vertex.outDeg() << " " << vertex.inDeg() << std::endl;
                    for (const Edge<htype> & e : vertex.getOutgoing()) {
                        std::cout << e.seq << std::endl;
                    }
                    VERIFY(false);
                }
                const Edge<htype> &edge = vertex.getOutgoing(seq[kmer.pos + hasher_.k]);
                Segment<Edge<htype>> seg(edge, 0, std::min(seq.size() - kmer.pos - hasher_.k, edge.size()));
                res.emplace_back(seg);
            }
        }
        return {prestart, res};
    }

    std::vector<std::pair<const Edge<htype> *, size_t>> carefulAlign(const Sequence & seq) const {
        std::vector<KWH<htype>> kmers = extractVertexPositions(seq);
        std::vector<std::pair<const Edge<htype>*, size_t>> res;
        if(kmers.size() == 0) {
            return res;
        }
        std::vector<const Vertex<htype> *> vertices;
        for(size_t i = 0; i < kmers.size(); i++) {
            vertices.push_back(&getVertex(kmers[i]));
        }
        for(size_t i = 0; i + 1 < kmers.size(); i++) {
            Sequence relevant= seq.Subseq(kmers[i].pos, kmers[i + 1].pos + hasher_.k);
            const Edge<htype> * best = nullptr;
            size_t best_score = 0;
            for(const Edge<htype> &edge : vertices[i]->getOutgoing()) {
                if (edge.end() == vertices[i + 1]) {
                    size_t score = std::min(
                        std::max(   hasher_.k,
                                    edge.common(relevant.Subseq(hasher_.k)) +
                                        vertices[i]->rcEdge(edge).common((!relevant).Subseq(hasher_.k))),
                        edge.size());
                    if(score > best_score) {
                        best_score = score;
                        best = &edge;
                    }
                }
            }
            if(best != nullptr) {
                res.emplace_back(best, best_score);
            }
        }
        return std::move(res);
    }

    void processRead(const Sequence & seq) {
        std::vector<KWH<htype>> kmers = extractVertexPositions(seq);
        VERIFY(kmers.size() > 0);
        std::vector<Vertex<htype> *> vertices;
        for(size_t i = 0; i < kmers.size(); i++) {
            vertices.emplace_back(&getVertex(kmers[i]));
            if(i == 0 || vertices[i] != vertices[i - 1]){
                vertices.back()->setSequence(kmers[i].getSeq());
                vertices.back()->incCoverage();
            }
        }
        for(size_t i = 0; i + 1 < vertices.size(); i++) {
//            TODO: if too memory heavy save only some of the labels
            VERIFY(kmers[i].pos + hasher_.k <= seq.size())
            if (i > 0 && vertices[i] == vertices[i - 1] && vertices[i] == vertices[i + 1] &&
                (kmers[i].pos - kmers[i - 1].pos == kmers[i + 1].pos - kmers[i].pos) &&
                kmers[i + 1].pos - kmers[i].pos < hasher_.k) {
                continue;
            }
            vertices[i]->addEdge(Edge<htype>(vertices[i + 1], Sequence(seq.Subseq(kmers[i].pos + hasher_.k,
                    kmers[i + 1].pos + hasher_.k).str())));
            vertices[i + 1]->rc().addEdge(Edge<htype>(&vertices[i]->rc(),
                    !Sequence(seq.Subseq(kmers[i].pos, kmers[i + 1].pos).str())));
        }
        if (kmers.front().pos > 0) {
            vertices.front()->rc().addEdge(Edge<htype>(nullptr, !(seq.Subseq(0, kmers[0].pos))));
        }
        if (kmers.back().pos + hasher_.k < seq.size()) {
            vertices.back()->addEdge(Edge<htype>(nullptr, seq.Subseq(kmers.back().pos + hasher_.k, seq.size())));
        }
    }

    size_t size() const {
        return v.size();
    }

    void printStats(logging::Logger &logger) const {
        std::vector<size_t> arr(10);
        size_t isolated = 0;
        size_t isolatedSize = 0;
        size_t n11 = 0;
        size_t n01 = 0;
        size_t e = 0;
        std::vector<size_t> inout(25);
        for(auto &val: v) {
            const Vertex<htype> &tmp = val.second;
            if (tmp.inDeg() == 0 && tmp.outDeg() == 1) {
                std::vector<Edge<htype>> path = tmp.walkForward(tmp.getOutgoing()[0]);
                if (path.back().end() != nullptr && path.back().end()->outDeg() == 0 && path.back().end()->inDeg() == 1) {
                    isolated += 1;
                    for(auto & edge : path) {
                        isolatedSize += edge.size();
                    }
                    isolatedSize += hasher().k;
                }
            }
            if (tmp.inDeg() == 1 && tmp.outDeg() == 0) {
                std::vector<Edge<htype>> path = tmp.rc().walkForward(tmp.rc().getOutgoing()[0]);
                if (path.back().end() != nullptr && path.back().end()->outDeg() == 0 && path.back().end()->inDeg() == 1) {
                    isolated += 1;
                    for(auto & edge : path) {
                        isolatedSize += edge.size();
                    }
                    isolatedSize += hasher().k;
                }
            }
            e == tmp.outDeg() + tmp.inDeg();
            arr[std::min(arr.size() - 1, tmp.outDeg())] += 1;
            arr[std::min(arr.size() - 1, tmp.inDeg())] += 1;
            inout[std::min<size_t>(4u, tmp.outDeg()) * 5 + std::min<size_t>(4u, tmp.inDeg())] += 1;
            inout[std::min<size_t>(4u, tmp.inDeg()) * 5 + std::min<size_t>(4u, tmp.outDeg())] += 1;
            if(tmp.outDeg() == 1 && tmp.inDeg() == 0) {
                Vertex<htype> & tmp1 = tmp.getOutgoing()[0].end()->rc();
                VERIFY(tmp.getOutgoing()[0].end() != nullptr);
            }
            if(tmp.outDeg() == 0 && tmp.inDeg() == 1) {
                Vertex<htype> &tmp1 = tmp.rc().getOutgoing()[0].end()->rc();
                VERIFY(tmp.rc().getOutgoing()[0].end()!= nullptr);
            }
            if (tmp.inDeg() == 1 && tmp.outDeg() == 1) {
                n11 += 1;
            }
            if (tmp.inDeg() + tmp.outDeg() == 1) {
                n01 += 1;
            }
            for(const Edge<htype> &edge : tmp.getOutgoing()) {
                e += 1;
            }
            for(const Edge<htype> &edge : tmp.rc().getOutgoing()) {
                e += 1;
            }
        }
        logger << "Graph statistics:" << std::endl;
        logger << "Total edges: " << e / 2 << std::endl;
        logger << "Total vertices: " << v.size() << std::endl;
        logger << "Number of end vertices: " << n01 << std::endl;
        logger << "Number of unbranching vertices: " << n11 << std::endl;
        logger << "Number of isolated edges " << isolated << " " << isolatedSize << std::endl;
        logger << "Distribution of degrees:" << std::endl;
        for(size_t i = 0; i < arr.size(); i++) {
            logger.noTimeSpace() << i << " " << arr[i] << std::endl;
        }
        logger << "Distribution of in/out degrees:" << std::endl;
        for(size_t i = 0; i < inout.size(); i++) {
            logger.noTimeSpace() << inout[i] << " ";
            if(i % 5 == 4)
                logger.noTimeSpace() << std::endl;
        }
    }

    void printCoverageStats(logging::Logger &logger) const {
        std::vector<size_t> cov(100);
        std::vector<size_t> cov_tips(100);
        std::vector<std::vector<size_t>> cov_ldist(100);
        std::vector<std::vector<size_t>> cov_ldist_tips(100);
        for(size_t i = 0; i < cov_ldist.size(); i++) {
            cov_ldist[i].resize(30);
        }
        for(size_t i = 0; i < cov_ldist.size(); i++) {
            cov_ldist_tips[i].resize(30);
        }
        std::vector<size_t> covLen(100);
        std::vector<size_t> covLen_tips(100);
        for(auto &val: v) {
            const Vertex<htype> &tmp = val.second;
            for(const Edge<htype> &edge : tmp.getOutgoing()) {
                cov[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                covLen[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
                cov_ldist[std::min(size_t(edge.getCoverage()), cov.size() - 1)][std::min(edge.size() / 50, cov_ldist[0].size() - 1)] += 1;
                if(edge.getTipSize() < hasher_.k * 2) {
                    cov_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                    covLen_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
                    cov_ldist_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)][std::min(edge.size() / 50, cov_ldist[0].size() - 1)] += 1;
                }
            }
            for(const Edge<htype> &edge : tmp.rc().getOutgoing()) {
                cov[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                covLen[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
                cov_ldist[std::min(size_t(edge.getCoverage()), cov.size() - 1)][std::min(edge.size() / 50, cov_ldist[0].size() - 1)] += 1;
                if(edge.getTipSize() < hasher_.k * 2) {
                    cov_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                    covLen_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
                    cov_ldist_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)][std::min(edge.size() / 50, cov_ldist[0].size() - 1)] += 1;
                }
            }
        }
        logger << "Distribution of coverages:" << std::endl;
        for(size_t i = 0; i < cov.size(); i++) {
            logger.noTimeSpace() << i << " " << cov[i] << " " << covLen[i];
            for(size_t val : cov_ldist[i]) {
                logger.noTimeSpace() << " " << val;
            }
            std::cout << std::endl;
            logger.noTimeSpace() << i << " " << cov_tips[i] << " " << covLen_tips[i];
            for(size_t val : cov_ldist_tips[i]) {
                logger.noTimeSpace() << " " << val;
            }
            std::cout << std::endl;
        }
    }

    typename std::unordered_map<htype, Vertex<htype>>::iterator begin() {
        return v.begin();
    }

    typename std::unordered_map<htype, Vertex<htype>>::iterator end() {
        return v.end();
    }

    typename std::unordered_map<htype, Vertex<htype>>::const_iterator begin() const {
        return v.begin();
    }

    typename std::unordered_map<htype, Vertex<htype>>::const_iterator end() const {
        return v.end();
    }

    void removeIsolated() {
        std::unordered_map<htype, Vertex<htype>> newv;
        for(auto & item : v) {
            if(item.second.outDeg() != 0 || item.second.inDeg() != 0) {
                newv.emplace(item.first, std::move(item.second));
            }
        }
        std::swap(v, newv);
    }

    void removeMarked() {
        std::unordered_map<htype, Vertex<htype>> newv;
        for(auto & item : v) {
            if(!item.second.marked()) {
                newv.emplace(item.first, std::move(item.second));
            }
        }
        std::swap(v, newv);
    }

    void printFasta(std::ostream &out) const {
        size_t cnt = 0;
        for(const auto &it : v) {
            const Vertex<htype> &vertex = it.second;
            VERIFY(!vertex.seq.empty());
            for(size_t i = 0; i < vertex.outDeg(); i++) {
                Sequence tmp = vertex.seq + vertex.getOutgoing()[i].seq;
                Vertex<htype> &end = *vertex.getOutgoing()[i].end();
                out << ">" << cnt << "_" << vertex.hash() << "1_" << end.hash() << int(end.isCanonical()) << std::endl;
                cnt++;
                out << tmp.str() << "\n";
            }
            const Vertex<htype> &rcvertex = vertex.rc();
            for(size_t i = 0; i < rcvertex.outDeg(); i++) {
                Sequence tmp = rcvertex.seq + rcvertex.getOutgoing()[i].seq;
                Vertex<htype> &end = *rcvertex.getOutgoing()[i].end();
                out << ">" << cnt << "_" << rcvertex.hash() << "1_" << end.hash() << int(end.isCanonical()) << std::endl;
                cnt++;
                out << tmp.str() << "\n";
            }
        }
    }

    template<class Iterator>
    void fillSparseDBGEdges(Iterator begin, Iterator end, logging::Logger &logger, size_t threads, const size_t min_read_size) {
        typedef typename Iterator::value_type ContigType;
        logger << "Starting to fill edges" << std::endl;
        std::function<void(ContigType &)> task = [this, min_read_size](ContigType & contig) {
            Sequence seq = contig.makeCompressedSequence();
            if(seq.size() >= min_read_size)
                processRead(seq);
        };
        processRecords(begin, end, logger, threads, task);
        logger << "Sparse graph edges filled." << std::endl;
    }

    static SparseDBG<htype> loadDBGFromFasta(const io::Library &lib, RollingHash<htype> & hasher, logging::Logger &logger, size_t threads) {
        logger << "Loading graph from fasta" << std::endl;
        io::SeqReader reader(lib);
        ParallelRecordCollector<Sequence> sequences(threads);
        ParallelRecordCollector<htype> vertices(threads);
        std::function<void(StringContig &)> collect_task = [&sequences, &vertices, hasher] (StringContig &contig){
            Sequence seq = contig.makeCompressedSequence();
            KWH<htype> start(hasher, seq, 0);
            KWH<htype> end(hasher, !seq, 0);
            vertices.add(start.hash());
            vertices.add(end.hash());
            sequences.add(seq);
        };
        processRecords(reader.begin(), reader.end(), logger, threads, collect_task);
        SparseDBG<htype> res(vertices.begin(), vertices.end(), hasher);
        reader.reset();
        res.fillSparseDBGEdges(reader.begin(), reader.end(), logger, threads, hasher.k + 1);
        logger << "Finished loading graph" << std::endl;
        return std::move(res);
    }
};

template<class htype>
class Component {
private:
    SparseDBG<htype> & graph;
    std::unordered_set<htype> v;
    struct EdgeRec {
        Vertex<htype> * start;
        Vertex<htype> * end;
        size_t size;
        size_t cov;
    };

    size_t outDeg(const Vertex<htype> &vert, size_t min_cov) const {
        size_t res = 0;
        for(const Edge<htype> &edge : vert.getOutgoing()) {
            if(edge.getCoverage() >= min_cov) {
                res += 1;
            }
        }
        return res;
    }

    bool isUnbranching(const Vertex<htype> &vert, size_t min_cov) const {
        return v.find(vert.hash()) != v.end() && outDeg(vert, min_cov) == 1 && outDeg(vert.rc(), min_cov) == 1;
    }

    Edge<htype> &getOut(Vertex<htype> &vert, size_t min_cov) {
        for(Edge<htype> &edge : vert.getOutgoing()) {
            if(edge.getCoverage() >= min_cov) {
                return edge;
            }
        }
        VERIFY(false);
    }

    Path<htype> unbranching(Vertex<htype> &vert, Edge<htype> &edge, size_t minCov) {
        std::vector<Edge<htype>> res;
        res.push_back(edge);
        Vertex<htype> *cur = edge.end();
        while (cur != &vert && isUnbranching(*cur, minCov)) {
            res.push_back(getOut(*cur, minCov));
            cur = res.back().end();
        }
        return Path<htype>(vert, res);
    }

public:
    template<class I>
    Component(SparseDBG<htype> &_graph, I begin, I end) : graph(_graph), v(begin, end) {
    }

    template<class I>
    static Component<htype> neighbourhood(SparseDBG<htype> &graph, I begin, I end, size_t radius, size_t min_coverage = 0) {
        std::unordered_set<htype> v;
        std::priority_queue<std::pair<size_t, htype>> queue;
        while(begin != end) {
            queue.emplace(0, *begin);
            ++begin;
        }
        while(!queue.empty()) {
            std::pair<size_t, htype> val = queue.top();
            queue.pop();
            if(v.find(val.second) != v.end())
                continue;
            v.insert(val.second);
            if(val.first > radius)
                continue;
            Vertex<htype> &vert = graph.getVertex(val.second);
            for(Edge<htype> & edge : vert.getOutgoing()) {
                if(edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
            for(Edge<htype> & edge : vert.rc().getOutgoing()) {
                if(edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
        }
        return Component<htype>(graph, v.begin(), v.end());
    }

    void printEdge(std::ostream &os, Vertex<htype> & start, Edge<htype> &edge) {
        Vertex<htype> &end = *edge.end();
        os << "\"";
        if (!start.isCanonical())
            os << "-";
        os << start.hash() % 10000 << "\" -> \"";
        if (!end.isCanonical())
            os << "-";
        os << end.hash() % 10000  << "\" [label=\"" << edge.size() << "(" << edge.getCoverage() << ")\"]\n";
    }

    void printEdge(std::ostream &os, Path<htype> & path) {
        size_t len = 0;
        size_t cov = 0;
        for(size_t i = 0; i < path.size(); i++) {
            len += path[i].size();
            cov += path[i].intCov();
        }
        Vertex<htype> &start = path.start();
        Vertex<htype> &end = *path.back().end();
        os << "\"";
        if (!start.isCanonical())
            os << "-";
        os << start.hash() % 10000 << "\" -> \"";
        if (!end.isCanonical())
            os << "-";
        os << end.hash() % 10000  << "\" [label=\"" << len << "(" << double(cov) / len << ")\"]\n";
    }

    size_t size() const {
        return v.size();
    }

    void printDot(std::ostream &os, size_t min_cov = 0) {
        os << "digraph {\nnodesep = 0.5;\n";
        for(htype vid : v) {
            Vertex<htype> &start = graph.getVertex(vid);
            for(Edge<htype> &edge : start.getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Vertex<htype> &end = *edge.end();
                printEdge(os, start, edge);
                if(v.find(end.hash()) == v.end()) {
                    printEdge(os, end.rc(), start.rcEdge(edge));
                }
            }
            for(Edge<htype> &edge : start.rc().getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Vertex<htype> &end = *edge.end();
                printEdge(os, start.rc(), edge);
                if(v.find(end.hash()) == v.end()) {
                    printEdge(os, end.rc(), start.rc().rcEdge(edge));
                }
            }
        }
        os << "}\n";
    }
    void printCompressedDot(std::ostream &os, size_t min_cov = 0) {
        os << "digraph {\nnodesep = 0.5;\n";
        for(htype vid : v) {
            Vertex<htype> &start = graph.getVertex(vid);
            if(isUnbranching(start, min_cov))
                continue;
            for(Edge<htype> &edge : start.getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Path<htype> path = unbranching(start, edge, min_cov);
                printEdge(os, path);
                Vertex<htype> &end = *path.back().end();
                if(v.find(end.hash()) == v.end()) {
                    Path<htype> rcpath = path.RC();
                    printEdge(os, rcpath);
                }
            }
            for(Edge<htype> &edge : start.rc().getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Path<htype> path = unbranching(start.rc(), edge, min_cov);
                printEdge(os, path);
                Vertex<htype> &end = *path.back().end();
                if(v.find(end.hash()) == v.end()) {
                    Path<htype> rcpath = path.RC();
                    printEdge(os, rcpath);
                }
            }
        }
        os << "}\n";
    }
};

template<typename htype, class Iterator>
void fillCoverage(SparseDBG<htype> &sdbg, logging::Logger &logger, Iterator begin, Iterator end, size_t threads,
                        const RollingHash<htype> &hasher, const size_t min_read_size) {
    typedef typename Iterator::value_type ContigType;
    logger << "Starting to fill edge coverages" << std::endl;
    ParallelRecordCollector<size_t> lens(threads);
    std::function<void(ContigType &)> task = [&sdbg, &lens, min_read_size](ContigType & contig) {
        Sequence seq = std::move(contig.makeCompressedSequence());
        if(seq.size() >= min_read_size) {
            GraphAlignment<htype> path = sdbg.align(seq);
            lens.add(path.size());
            for(Segment<Edge<htype>> &seg : path) {
                seg.contig().incCov(seg.size());
            }
            path = sdbg.align(!seq);
            for(Segment<Edge<htype>> &seg : path) {
                seg.contig().incCov(seg.size());
            }
        }
    };
    processRecords(begin, end, logger, threads, task);
    logger << "Edge coverage calculated." << std::endl;
    std::vector<size_t> lens_distr(1000);
    for(size_t l : lens) {
        lens_distr[std::min(l, lens_distr.size() - 1)] += 1;
    }
    logger << "Distribution of path sizes." << std::endl;
    for(size_t i = 0; i < lens_distr.size(); i++)
        std::cout << i << " " << lens_distr[i] << std::endl;
}

template<typename htype>
SparseDBG<htype> constructSparseDBGFromReads(logging::Logger & logger, const io::Library &reads_file, size_t threads, const RollingHash<htype> &hasher,
                                       const std::vector<htype> &hash_list, const size_t w) {
    logger << "Starting construction of sparse de Bruijn graph" << std::endl;
    SparseDBG<htype> sdbg(hash_list.begin(), hash_list.end(), hasher);
    logger << "Vertex map constructed." << std::endl;
    io::SeqReader reader(reads_file);
    sdbg.fillSparseDBGEdges(reader.begin(), reader.end(), logger, threads, w + hasher.k - 1);
    return std::move(sdbg);
}


template<typename htype>
void tieTips(logging::Logger &logger, SparseDBG<htype> &sdbg, size_t w, size_t threads) {
    logger << " Collecting tips " << std::endl;
//    TODO reduce memory consumption!! A lot of duplicated k-mer storing
    ParallelRecordCollector<Sequence> old_edges(threads);
    ParallelRecordCollector<htype> new_minimizers(threads);
    std::function<void(std::pair<const htype, Vertex<htype>> &)> task =
            [&sdbg, &old_edges, &new_minimizers](std::pair<const htype, Vertex<htype>> & pair) {
        Vertex<htype> &rec = pair.second;
        VERIFY(!rec.seq.empty());
        for (size_t i = 0; i < rec.getOutgoing().size(); i++) {
            const auto &ext = rec.getOutgoing()[i];
            Sequence seq = rec.seq + ext.seq;
            old_edges.add(seq);
            if (ext.end() == nullptr) {
                KWH<htype> kwh(sdbg.hasher(), seq, ext.size());
                new_minimizers.emplace_back(kwh.hash());
            }
        }
        const Vertex<htype> &rec1 = rec.rc();
        for (size_t i = 0; i < rec1.getOutgoing().size(); i++) {
            const auto &ext = rec1.getOutgoing()[i];
            Sequence seq = rec1.seq + ext.seq;
            old_edges.add(seq);
            if (ext.end() == nullptr) {
                KWH<htype> kwh(sdbg.hasher(), seq, ext.size());
                new_minimizers.emplace_back(kwh.hash());
            }
        }
        rec.clear();
    };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);

//#pragma omp parallel default(none) shared(sdbg, old_edges, new_minimizers, logger)
//    {
//#pragma omp single
//        {
//            for (auto &it: sdbg) {
//                Vertex<htype> &rec = it.second;
//                VERIFY(!rec.seq.empty());
//#pragma omp task default(none) shared(sdbg, old_edges, new_minimizers, rec, logger)
//                {
//                    for (size_t i = 0; i < rec.getOutgoing().size(); i++) {
//                        const auto &ext = rec.getOutgoing()[i];
//                        Sequence seq = rec.seq + ext.seq;
//                        old_edges.add(seq);
//                        if (ext.end() == nullptr) {
//                            KWH<htype> kwh(sdbg.hasher(), seq, ext.size());
//                            new_minimizers.emplace_back(kwh.hash());
//                        }
//                    }
//                    const Vertex<htype> &rec1 = rec.rc();
//                    for (size_t i = 0; i < rec1.getOutgoing().size(); i++) {
//                        const auto &ext = rec1.getOutgoing()[i];
//                        Sequence seq = rec1.seq + ext.seq;
//                        old_edges.add(seq);
//                        if (ext.end() == nullptr) {
//                            KWH<htype> kwh(sdbg.hasher(), seq, ext.size());
//                            new_minimizers.emplace_back(kwh.hash());
//                        }
//                    }
//                    rec.clear();
//                }
//            }
//        }
//    }
    logger << "Added " << new_minimizers.size() << " artificial minimizers from tips." << std::endl;
    logger << "Collected " << old_edges.size() << " old edges." << std::endl;
    for(auto it = new_minimizers.begin(); it != new_minimizers.end(); ++it) {
        sdbg.addVertex(*it);
    }
    logger << "New minimizers added to sparse graph." << std::endl;
    logger << "Refilling graph with edges." << std::endl;
    sdbg.fillSparseDBGEdges(old_edges.begin(), old_edges.end(), logger, threads, sdbg.hasher().k + 1);
    logger << "Finished fixing sparse de Bruijn graph." << std::endl;
}

template<typename htype>
void UpdateVertexTips(Vertex<htype> &rec, ParallelRecordCollector<Vertex<htype> *> &queue) {
    bool ok = true;
    for (const Edge<htype> &edge : rec.getOutgoing()) {
        if (edge.getTipSize() == size_t(-1)) {
            edge.updateTipSize();
        }
        if (edge.getTipSize() == size_t(-1)) {
            ok = false;
        }
    }
    if(ok && rec.inDeg() == 1) {
        queue.add(&(rec.rc().getOutgoing()[0].end()->rc()));
    }
}

template<typename htype>
void findTips(logging::Logger &logger, SparseDBG<htype> &sdbg, size_t threads) {
    logger << " Finding tips " << std::endl;
//    TODO reduce memory consumption!! A lot of duplicated k-mer storing
    ParallelRecordCollector<Vertex<htype> *> queue(threads);
#pragma omp parallel default(none) shared(sdbg, logger, queue)
    {
#pragma omp single
        {
            for (auto &it: sdbg) {
                Vertex<htype> &rec = it.second;
                VERIFY(!rec.seq.empty());
#pragma omp task default(none) shared(sdbg, rec, logger, queue)
                {
                    UpdateVertexTips(rec, queue);
                    UpdateVertexTips(rec.rc(), queue);
                }
            }
        }
    }
    logger << "Found initial tips. Looking for iterative tips" << std::endl;
    size_t cnt = 0;
    while(!queue.empty()) {
        logger << "Iteration " << cnt << ". Queue size " << queue.size() << std::endl;
        std::vector<Vertex<htype> *> prev_queue = queue.collectUnique();
        queue.clear();
#pragma omp parallel default(none) shared(sdbg, logger, prev_queue, queue)
        {
#pragma omp single
            {
                for (auto &it: prev_queue) {
                    Vertex<htype> &rec = *it;
                    VERIFY(!rec.seq.empty());
#pragma omp task default(none) shared(sdbg, rec, logger, queue)
                    {
                        UpdateVertexTips(rec, queue);
                    }
                }
            }
        }
    }
    logger << "Tip finding finished" << std::endl;
}


template<typename htype>
void mergeLoop(Vertex<htype> &start, std::vector<Edge<htype>> &path) {
    if(path.size() % 2 == 0 && *path[path.size() / 2].end() == start.rc()) {
        path =std::vector<Edge<htype>>(path.begin(), path.begin() + path.size() / 2);
    }
    start.rc().lock();
    for(const Edge<htype> &e : path) {
        e.end()->lock();
        e.end()->rc().lock();
    }
    Sequence newSeq(start.pathSeq(path));
    for(const Edge<htype> &e : path) {
        e.end()->mark();
    }
    start.addEdgeLockFree(Edge<htype>(path.back().end(), newSeq.Subseq(start.seq.size())));
    path.back().end()->rc().addEdgeLockFree(Edge<htype>(&start.rc(), (!newSeq).Subseq(start.seq.size())));
    for(const Edge<htype> &e : path) {
        e.end()->unlock();
        e.end()->rc().unlock();
    }
}

template<class htype>
void MergeEdge(SparseDBG<htype> &sdbg, Vertex<htype> &start, const Edge<htype> &edge) {
    std::vector<Edge<htype>> path = start.walkForward(edge);
    Vertex<htype> &end = path.back().end()->rc();
    if (path.size() > 1 && end.hash() >= start.hash()) {
        VERIFY(start.seq.size() > 0)
        VERIFY(end.seq.size() > 0);
        Sequence newSeq(start.pathSeq(path));
        if (start != end)
            end.lock();
        start.addEdgeLockFree(Edge<htype>(&end.rc(), newSeq.Subseq(start.seq.size())));
        end.addEdgeLockFree(Edge<htype>(&start.rc(), (!newSeq).Subseq(start.seq.size())));
        for(size_t i = 0; i + 1 < path.size(); i++) {
//            path[i].end()->clear();
            path[i].end()->mark();
            path[i].end()->rc().mark();
        }
        if (start != end)
            end.unlock();
    }
}

template<class htype>
void mergeLinearPaths(logging::Logger & logger, SparseDBG<htype> &sdbg, size_t threads) {
    logger << "Merging linear unbranching paths" << std::endl;
    std::function<void(std::pair<const htype, Vertex<htype>> &)> task =
            [&sdbg](std::pair<const htype, Vertex<htype>> & pair) {
                Vertex<htype> &start = pair.second;
                if (!start.isJunction())
                    return;
                start.lock();
                for (const Edge<htype> &edge: start.getOutgoing()) {
                    MergeEdge(sdbg, start, edge);
                }
                start.unlock();
                start.rc().lock();
                for (const Edge<htype> &edge: start.rc().getOutgoing()) {
                    MergeEdge(sdbg, start.rc(), edge);
                }
                start.rc().unlock();
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);

//#pragma omp parallel default(none) shared(sdbg)
//    {
//#pragma omp single
//        {
//            for (auto &it: sdbg) {
//                Vertex<htype> &start = it.second;
//                if (!start.isJunction())
//                    continue;
//#pragma omp task default(none) shared(start, sdbg)
//                {
//                    start.lock();
//                    for (const Edge<htype> &edge: start.getOutgoing()) {
//                        MergeEdge(sdbg, start, edge);
//                    }
//                    start.unlock();
//                    start.rc().lock();
//                    for (const Edge<htype> &edge: start.rc().getOutgoing()) {
//                        MergeEdge(sdbg, start.rc(), edge);
//                    }
//                    start.rc().unlock();
//                }
//            }
//        }
//    }
    logger << "Finished merging linear unbranching paths" << std::endl;
}

template<class htype>
void mergeCyclicPaths(logging::Logger & logger, SparseDBG<htype> &sdbg, size_t threads) {
    logger << "Merging cyclic paths" << std::endl;
    std::function<void(std::pair<const htype, Vertex<htype>> &)> task =
            [&sdbg](std::pair<const htype, Vertex<htype>> & pair) {
                Vertex<htype> &start = pair.second;
                start.lock();
                if(start.isJunction() || start.marked()) {
                    start.unlock();
                    return;
                }
                std::vector<Edge<htype>> path = start.walkForward(start.getOutgoing()[0]);
                if (*path.back().end() == start) {
                    bool ismin = true;
                    for (const Edge<htype> &e : path) {
                        if (e.end()->hash() < start.hash()) {
                            ismin = false;
                            break;
                        }
                    }
                    if(ismin) {
                        mergeLoop(start, path);
                        start.unlock();
                    }
                }
                start.unlock();
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);

//#pragma omp parallel default(none) shared(sdbg)
//    {
//#pragma omp single
//        {
//            for( auto &it: sdbg) {
//                Vertex<htype> &start =  it.second;
//#pragma omp task default(none) shared(start, sdbg)
//                {
//                    start.lock();
//                    if(start.inDeg() == 1 && start.outDeg() == 1) {
//                        std::vector<Edge<htype>> path = start.walkForward(start.getOutgoing()[0]);
//                        if (*path.back().end() == start) {
//                            bool ismin = true;
//                            for (const Edge<htype> &e : path) {
//                                if (e.end()->hash() < start.hash()) {
//                                    ismin = false;
//                                    break;
//                                }
//                            }
//                            if(ismin) {
//                                start.unlock();
//                                mergeLoop(start, path);
//                            }
//                        }
//                    }
//                    start.unlock();
//                }
//            }
//        }
//    }
    logger << "Finished merging cyclic paths" << std::endl;
}

template<class htype>
void mergeAll(logging::Logger & logger, SparseDBG<htype> &sdbg, size_t threads) {
    mergeLinearPaths(logger, sdbg, threads);
//    sdbg.checkConsistency(threads, logger);
    mergeCyclicPaths(logger, sdbg, threads);
//    sdbg.checkConsistency(threads, logger);
    logger << "Removing isolated vertices" << std::endl;
    sdbg.removeMarked();
    logger << "Finished removing isolated vertices" << std::endl;
    sdbg.checkConsistency(threads, logger);
}