//
// Created by anton on 7/22/20.
//

#pragma once
#include "sequences/sequence.hpp"
#include "sequences/seqio.hpp"
#include "common/omp_utils.hpp"
#include "common/logging.hpp"
#include "rolling_hash.hpp"
#include <vector>
#include <numeric>
#include <unordered_map>

template<class htype>
class Vertex;

template<class htype>
class Edge {
private:
    Sequence seq_;
    Vertex<htype> *end_;
    mutable size_t cov;
public:
    friend class Vertex<htype>;

    Edge(Vertex<htype> *_end, const Sequence &_seq) :
            seq_(_seq), end_(_end), cov(0) {
    }

    Vertex<htype> *end() const {
        return end_;
    }

    size_t size() const {
        return seq_.size();
    }

    double getCoverage() const {
        return double(cov) / size();
    }

    void incCov(size_t val) const {
#pragma omp atomic
        cov += val;
    }

    const Sequence & seq() const {
        return seq_;
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

public:
    Sequence seq;

    void clear() {
        outgoing_.clear();
        rc_->outgoing_.clear();
    }

    explicit Vertex(htype hash = 0, Vertex *_rc = nullptr) : hash_(hash), rc_(_rc) {
        if (rc_ == nullptr) {
            rc_ = new Vertex<htype>(hash, this);
        }
        omp_init_lock(&writelock);
    }

    Vertex(Vertex &&other) noexcept : rc_(other.rc_), hash_(other.hash_) {
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

    Vertex(const Vertex &) = delete;

    ~Vertex() {
        if(rc_ != nullptr) {
            rc_->rc_ = nullptr;
            free(rc_);
        }
        rc_ = nullptr;
    }

    void checkConsistency() const {
        for(const Edge<htype> & edge : outgoing_) {
            if(edge.end() != nullptr) {
                if(this->rcEdge(edge).end() != &(this->rc())) {
                    std::cout << this << " " << seq << " " << edge.seq() << " " << rcEdge(edge).end() << " " << &(this->rc()) <<std::endl;
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

    Edge<htype>& rcEdge(const Edge<htype> & edge) {
        Vertex &vend = edge.end()->rc();
        char c;
        if(edge.size() > seq.size()) {
            c = (!edge.seq())[seq.size()];
        } else {
            c = (!seq)[seq.size() - edge.size()];
        }
        return vend.getOutgoing(c);
    }

    const Edge<htype>& rcEdge(const Edge<htype> & edge) const {
        const Vertex &vend = edge.end()->rc();
        char c;
        if(edge.size() > seq.size()) {
            c = (!edge.seq())[seq.size()];
        } else {
            c = (!seq)[seq.size() - edge.size()];
        }
        return vend.getOutgoing(c);
    }

    Sequence pathSeq(const std::vector<Edge<htype>> & path) const {
        SequenceBuilder sb;
        sb.append(seq);
        for(const Edge<htype> &e : path) {
            sb.append(e.seq());
        }
        return sb.BuildSequence();
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
            unlock();
        }
//        else {
//            VERIFY(_seq == seq);
//        }
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
                if (edge.seq() == e.seq().Subseq(0, edge.size())) {
                    return;
                }
            } else if (edge.seq().Subseq(0, e.size()) == e.seq()) {
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
            if(edge.seq()[0] == c) {
                return edge;
            }
        }
        VERIFY(false);
        return getOutgoing()[0];
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
};

template<typename htype>
class SparseDBG {
private:
//    TODO: replace with perfect hash map? It is parallel, maybe faster and compact.
    std::unordered_map<htype, Vertex<htype>> v;
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

public:

    SparseDBG(const std::vector<htype> &_hashs, RollingHash<htype> _hasher) : hasher_(_hasher) {
        std::vector<htype> hashs(_hashs);
        for(htype hash: hashs) {
            addVertex(hash);
        }
    }

    SparseDBG(SparseDBG<htype> &&other) noexcept = default;

    void checkConsistency() {
        std::cout << "Checking consistency" << std::endl;
        for(const auto & it : v) {
            const Vertex<htype> & vert = it.second;
            vert.checkConsistency();
            vert.rc().checkConsistency();
        }
        std::cout << "Consistency check success" << std::endl;
    }

    const RollingHash<htype> &hasher() const {
        return hasher_;
    }

    void addVertex(htype h) {
        v.emplace(h, Vertex<htype>(h));
        VERIFY(v.find(h) != v.end());
    }

    Vertex<htype> &getVertex(const KWH<htype> &kwh) {
        VERIFY(v.find(kwh.hash()) != v.end());
        if(kwh.isCanonical()) {
            return v[kwh.hash()];
        } else {
            return v[kwh.hash()].rc();
        }
    }

    const Vertex<htype> &getVertex(const KWH<htype> &kwh) const {
        VERIFY(v.find(kwh.hash()) != v.end());
        if(kwh.isCanonical()) {
            return v.find(kwh.hash())->second;
        } else {
            return v.find(kwh.hash())->second.rc();
        }
    }

    std::vector<Segment<Edge<htype>>> align(const Sequence & seq) const {
        std::vector<KWH<htype>> kmers = extractVertexPositions(seq);
        std::vector<Segment<Edge<htype>>> res;
        if(kmers.size() == 0) {
            return res;
        }
        if (kmers.front().pos > 0) {
            const Vertex<htype> &rcstart = getVertex(kmers.front()).rc();
            const Edge<htype> &edge = rcstart.rcEdge(rcstart.getOutgoing(seq[kmers.front().pos - 1] ^ 3));
            res.emplace_back(edge, edge.size() - kmers.front().pos, edge.size());
        }
        for(const KWH<htype> & kmer : kmers) {
            if (kmer.pos + hasher_.k < seq.size()) {
                const Vertex<htype> &vertex = getVertex(kmer);
                const Edge<htype> &edge = vertex.getOutgoing(seq[kmer.pos + hasher_.k]);
                res.emplace_back(edge, 0, std::min(seq.size() - kmer.pos, edge.size()));
            }
        }
        return std::move(res);
    }

//    Add all edges (including hanging) from read. All Sequences are copied and read should be discarded after.
    void processRead(const Sequence & seq) {
        std::vector<KWH<htype>> kmers = extractVertexPositions(seq);
        VERIFY(kmers.size() > 0);
        std::vector<Vertex<htype> *> vertices;
        for(size_t i = 0; i < kmers.size(); i++) {
            vertices.emplace_back(&getVertex(kmers[i]));
            if(i == 0 || vertices[i] != vertices[i - 1]){
                vertices.back()->setSequence(kmers[i].getSeq());
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

//    Add proper edges from disjointig. All sequences are subsequences of disjointig and thus not copied.
    void processDisjointig(const Sequence & seq) {
        std::vector<KWH<htype>> kmers = extractVertexPositions(seq);
        std::vector<Vertex<htype> *> vertices;
        for(size_t i = 0; i < kmers.size(); i++) {
            vertices.emplace_back(&getVertex(kmers[i]));
            if(i == 0 || vertices[i] != vertices[i - 1])
                vertices.back()->setSequence(kmers[i].getSeq());
        }
        for(size_t i = 0; i + 1 < vertices.size(); i++) {
//            TODO: if too memory heavy save only some of the labels
            if (i > 0 && vertices[i] == vertices[i - 1] && vertices[i] == vertices[i + 1] &&
                (kmers[i].pos - kmers[i - 1].pos == kmers[i + 1].pos - kmers[i].pos) &&
                kmers[i + 1].pos - kmers[i].pos < hasher_.k) {
                continue;
            }
            vertices[i]->addEdge(Edge<htype>(vertices[i + 1], seq.Subseq(kmers[i].pos + hasher_.k, kmers[i + 1].pos + hasher_.k)));
            vertices[i + 1]->rc().addEdge(Edge<htype>(&vertices[i]->rc(), !(seq.Subseq(kmers[i].pos, kmers[i + 1].pos))));
        }
    }

    std::vector<Edge<htype>> walkForward(const Vertex<htype> &rec, const Edge<htype> & edge) const {
        std::vector<Edge<htype>> res;
        res.push_back(edge);
        Vertex<htype> *next = edge.end();
        while(next != nullptr && next != &rec && !next->isJunction()) {
            res.push_back(next->getOutgoing()[0]);
            next = next->getOutgoing()[0].end();
        }
        return res;
    }

    size_t size() const {
        return v.size();
    }

    void printStats(logging::Logger &logger) const {
        std::vector<size_t> arr(10);
        std::vector<size_t> cov(20);
        std::vector<size_t> covLen(20);
        size_t isolated = 0;
        size_t isolatedSize = 0;
        size_t n11 = 0;
        size_t n01 = 0;
        size_t e = 0;
        std::vector<size_t> inout(25);
        for(auto &val: v) {
            const Vertex<htype> &tmp = val.second;
            e == tmp.outDeg() + tmp.inDeg();
            arr[std::min(arr.size() - 1, tmp.outDeg())] += 1;
            arr[std::min(arr.size() - 1, tmp.inDeg())] += 1;
            inout[std::min<size_t>(4u, tmp.outDeg()) * 5 + std::min<size_t>(4u, tmp.inDeg())] += 1;
            inout[std::min<size_t>(4u, tmp.inDeg()) * 5 + std::min<size_t>(4u, tmp.outDeg())] += 1;
            if(tmp.outDeg() == 1 && tmp.inDeg() == 0) {
                Vertex<htype> & tmp1 = tmp.getOutgoing()[0].end()->rc();
                VERIFY(tmp.getOutgoing()[0].end() != nullptr);
                if (tmp1.outDeg() == 1 && tmp.inDeg() == 0) {
                    isolated += 1;
                    isolatedSize += tmp1.getOutgoing()[0].size();
                }
            }
            if(tmp.outDeg() == 0 && tmp.inDeg() == 1) {
                Vertex<htype> &tmp1 = tmp.rc().getOutgoing()[0].end()->rc();
                VERIFY(tmp.rc().getOutgoing()[0].end()!= nullptr);
                if (tmp1.outDeg() == 1 && tmp.inDeg() == 0) {
                    isolated += 1;
                    isolatedSize += tmp1.getOutgoing()[0].size();
                }
            }
            if (tmp.inDeg() == 1 && tmp.outDeg() == 1) {
                n11 += 1;
            }
            if (tmp.inDeg() + tmp.outDeg() == 1) {
                n01 += 1;
            }
            for(const Edge<htype> &edge : tmp.getOutgoing()) {
                e += 1;
                cov[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                covLen[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
            }
            for(const Edge<htype> &edge : tmp.rc().getOutgoing()) {
                e += 1;
                cov[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                covLen[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
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
        logger << "Distribution of coverages:" << std::endl;
        for(size_t i = 0; i < cov.size(); i++) {
            logger.noTimeSpace() << i << " " << cov[i] << " " << covLen[i] << std::endl;
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

    void printFasta(std::ostream &out) const {
        size_t cnt = 0;
        for(const auto &it : v) {
            const Vertex<htype> &vertex = it.second;
            for(size_t i = 0; i < vertex.outDeg(); i++) {
                Sequence tmp = vertex.seq + vertex.getOutgoing()[i].seq();
                if (tmp <= !tmp) {
                    out << ">" << cnt << "\n";
                    cnt++;
                    out << tmp.str() << "\n";
                }
            }
            for(size_t i = 0; i < vertex.inDeg(); i++) {
                Sequence tmp = vertex.rc().seq + vertex.rc().getOutgoing()[i].seq();
                if (tmp <= !tmp) {
                    out << ">" << cnt << "\n";
                    cnt++;
                    out << tmp.str() << "\n";
                }
            }
        }
    }
};

template<typename htype, class Iterator>
void fillSparseDBGEdges(SparseDBG<htype> &sdbg, logging::Logger &logger, Iterator begin, Iterator end, size_t threads,
                        const RollingHash<htype> &hasher, const size_t min_read_size) {
    logger << "Starting to fill edges" << std::endl;
    const size_t buffer_size = 1000000000;
    while(begin != end) {
        size_t tlen = 0;
        logger << "Starting new round" << std::endl;
        std::vector<Sequence> reads;
        reads.reserve(1000000);
#pragma omp parallel default(none) shared(hasher, begin, end, tlen, buffer_size, std::cout, logger, reads, sdbg, min_read_size)
        {
#pragma omp single
            {
                while (begin != end && tlen < buffer_size && reads.size() < 1000000) {
                    reads.push_back(*begin);
                    ++begin;
                    tlen += reads.back().size();
                    if(reads.back().size() < min_read_size) {
                        continue;
                    }
                    Sequence &read = reads.back();
#pragma omp task default(none) shared(read, sdbg, std::cout)
                    {
                        sdbg.processRead(read);
                    }
                }
                logger << tlen  << " nucleotides in " << reads.size() <<
                          " sequences were collected. Processing in progress  " << std::endl;
            }
        }
        reads.clear();
    }
    logger << "Sparse graph edges filled." << std::endl;
}

template<typename htype, class Iterator>
void fillCoverage(SparseDBG<htype> &sdbg, logging::Logger &logger, Iterator begin, Iterator end, size_t threads,
                        const RollingHash<htype> &hasher, const size_t min_read_size) {
    logger << "Starting to fill edge coverages" << std::endl;
    const size_t buffer_size = 1000000000;
    while(begin != end) {
        size_t tlen = 0;
        logger << "Starting new round" << std::endl;
        std::vector<Sequence> reads;
        reads.reserve(1000000);
#pragma omp parallel default(none) shared(hasher, begin, end, tlen, buffer_size, std::cout, logger, reads, sdbg, min_read_size)
        {
#pragma omp single
            {
                while (begin != end && tlen < buffer_size && reads.size() < 1000000) {
                    reads.push_back(*begin);
                    ++begin;
                    tlen += reads.back().size();
                    if(reads.back().size() < min_read_size) {
                        continue;
                    }
                    size_t index = reads.size() - 1;
#pragma omp task default(none) shared(reads, index, sdbg, std::cout)
                    {
//                        TODO make better for rc

                        std::vector<Segment<Edge<htype>>> path = sdbg.align(reads[index]);
                        for(Segment<Edge<htype>> &seg : path) {
                            seg.contig.incCov(seg.size());
                        }
                        path = sdbg.align(!reads[index]);
                        for(Segment<Edge<htype>> &seg : path) {
                            seg.contig.incCov(seg.size());
                        }
                    }
                }
                logger << tlen  << " nucleotides in " << reads.size() <<
                       " sequences were collected. Processing in progress  " << std::endl;
            }
        }
        reads.clear();
    }
    logger << "Sparse graph edges filled." << std::endl;
}

template<typename htype>
SparseDBG<htype> constructSparseDBGFromReads(logging::Logger & logger, const std::string &reads_file, size_t threads, const RollingHash<htype> &hasher,
                                       std::vector<htype> &&hash_list, const size_t w) {
    logger << "Starting construction of sparse de Bruijn graph" << std::endl;
    SparseDBG<htype> sdbg(std::move(hash_list), hasher);
    logger << "Vertex map constructed." << std::endl;
    io::SeqReader reader(io::SeqReader::CompressingReader(reads_file));
    fillSparseDBGEdges(sdbg, logger, reader.seqbegin(), reader.seqend(), threads, hasher, w + hasher.k - 1);
    return std::move(sdbg);
}


template<typename htype>
void tieTips(logging::Logger &logger, SparseDBG<htype> &sdbg, size_t w, size_t threads) {
    logger << " Collecting tips " << std::endl;
    ParallelRecordCollector<htype> new_minimizers(threads);
//    TODO reduce memory consumption!! A lot of duplicated k-mer storing
    ParallelRecordCollector<Sequence> old_edges(threads);
#pragma omp parallel default(none) shared(sdbg, old_edges, new_minimizers, logger)
    {
#pragma omp single
        {
            for (auto &it: sdbg) {
                Vertex<htype> &rec = it.second;
                VERIFY(!rec.seq.empty());
#pragma omp task default(none) shared(sdbg, old_edges, new_minimizers, rec, logger)
                {
                    for (size_t i = 0; i < rec.getOutgoing().size(); i++) {
                        const auto &ext = rec.getOutgoing()[i];
                        Sequence seq = rec.seq + ext.seq();
                        old_edges.add(seq);
                        if (ext.end() == nullptr) {
                            KWH<htype> kwh(sdbg.hasher(), seq, ext.seq().size());
                            new_minimizers.emplace_back(kwh.hash());
                        }
                    }
                    const Vertex<htype> &rec1 = rec.rc();
                    for (size_t i = 0; i < rec1.getOutgoing().size(); i++) {
                        const auto &ext = rec1.getOutgoing()[i];
                        Sequence seq = rec1.seq + ext.seq();
                        old_edges.add(seq);
                        if (ext.end() == nullptr) {
                            KWH<htype> kwh(sdbg.hasher(), seq, ext.seq().size());
                            new_minimizers.emplace_back(kwh.hash());
                        }
                    }
                    rec.clear();
                }
            }
        }
    }
    logger << "Added " << new_minimizers.size() << " artificial minimizers from tips." << std::endl;
    logger << "Collected " << old_edges.size() << " old edges." << std::endl;
    for(auto it = new_minimizers.begin(); it != new_minimizers.end(); ++it) {
        sdbg.addVertex(*it);
    }
    logger << "New minimizers added to sparse graph." << std::endl;
    logger << "Refilling graph with edges." << std::endl;
    fillSparseDBGEdges(sdbg, logger, old_edges.begin(), old_edges.end(), threads, sdbg.hasher(), sdbg.hasher().k + 1);
    logger << "Finished fixing sparse de Bruijn graph." << std::endl;
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
        e.end()->clear();
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
    std::vector<Edge<htype>> path = sdbg.walkForward(start, edge);
    Vertex<htype> &end = path.back().end()->rc();
    if (path.size() > 1 && end.hash() >= start.hash()) {
        VERIFY(start.seq.size() > 0)
        VERIFY(end.seq.size() > 0);
        Sequence newSeq(start.pathSeq(path));
        if (start != end)
            end.lock();
        start.addEdgeLockFree(Edge<htype>(&end.rc(), newSeq.Subseq(start.seq.size())));
        end.addEdgeLockFree(Edge<htype>(&start.rc(), (!newSeq).Subseq(start.seq.size())));
        if (start != end)
            end.unlock();
        for(size_t i = 0; i + 1 < path.size(); i++) {
            path[i].end()->clear();
        }
    }

}
template<class htype>
void mergeLinearPaths(logging::Logger & logger, SparseDBG<htype> &sdbg) {
    logger << "Merging linear unbranching paths" << std::endl;
#pragma omp parallel default(none) shared(sdbg)
    {
#pragma omp single
        {
            for (auto &it: sdbg) {
                Vertex<htype> &start = it.second;
                if (!start.isJunction())
                    continue;
#pragma omp task default(none) shared(start, sdbg)
                {
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
                }
            }
        }
    }
    logger << "Finished merging linear unbranching paths" << std::endl;
}

template<class htype>
void mergeCyclicPaths(logging::Logger & logger, SparseDBG<htype> &sdbg) {
    logger << "Merging cyclic paths" << std::endl;
#pragma omp parallel default(none) shared(sdbg)
    {
#pragma omp single
        {
            for( auto &it: sdbg) {
                Vertex<htype> &start =  it.second;
#pragma omp task default(none) shared(start, sdbg)
                {
                    start.lock();
                    if(start.inDeg() == 1 && start.outDeg() == 1) {
                        std::vector<Edge<htype>> path = sdbg.walkForward(start, start.getOutgoing()[0]);
                        if (*path.back().end() == start) {
                            bool ismin = true;
                            for (const Edge<htype> &e : path) {
                                if (e.end()->hash() < start.hash()) {
                                    ismin = false;
                                    break;
                                }
                            }
                            if(ismin) {
                                start.unlock();
                                mergeLoop(start, path);
                            }
                        }
                    }
                    start.unlock();
                }
            }
        }
    }
    logger << "Finished merging cyclic paths" << std::endl;
}

template<class htype>
void mergeAll(logging::Logger & logger, SparseDBG<htype> &sdbg) {
    mergeLinearPaths(logger, sdbg);
    sdbg.checkConsistency();
    mergeCyclicPaths(logger, sdbg);
    sdbg.checkConsistency();
    logger << "Removing isolated vertices" << std::endl;
    sdbg.removeIsolated();
    logger << "Finished removing isolated vertices" << std::endl;
    sdbg.checkConsistency();
}