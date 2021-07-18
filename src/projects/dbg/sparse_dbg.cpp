#include "sparse_dbg.hpp"
using namespace dbg;

size_t Edge::updateTipSize() const {
    size_t new_val = 0;
    if(extraInfo == size_t(-1) && end_->inDeg() == 1) {
        for (const Edge & other : *end_) {
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

Edge &Edge::rc() const {
    Vertex &vend = end_->rc();
    char c;
    size_t k = vend.seq.size();
    if(size() > k) {
        c = (!seq)[k];
    } else {
        c = (!start_->seq)[k - size()];
    }
    return vend.getOutgoing(c);
}

Edge &Edge::sparseRcEdge() const {
    Vertex &vend = end()->rc();
    VERIFY(start_->seq.size() > 0);
    for(Edge &candidate : vend) {
        if(*candidate.end() == start_->rc() && candidate.size() == size() &&
           (size() <= start_->seq.size() || candidate.seq.startsWith((!seq).Subseq(start_->seq.size())))) {
            return candidate;
        }
    }
    std::cout << start_->seq + seq << std::endl;
    std::cout << seq << std::endl;
    std::cout << vend.seq << std::endl;
    for(Edge &candidate : vend) {
        if(*candidate.end() == start_->rc() && candidate.size() == size() &&
           (size() <= seq.size() || candidate.seq.startsWith((!seq).Subseq(start_->seq.size())))) {
            std::cout << vend.seq + candidate.seq << std::endl;
        }
    }
    VERIFY(false);
    return vend[0];
}

void Edge::bindTip(Vertex &start, Vertex &end) {
    VERIFY(end_ == nullptr);
    end_ = &end;
    Sequence rcseq = !(start.seq + seq);
    end.rc().addEdgeLockFree(Edge(&end.rc(), &start.rc(), rcseq.Subseq(start.seq.size())));
}

size_t Edge::getTipSize() const {
    return extraInfo;
}

Vertex *Edge::start() const {
    return start_;
}

Vertex *Edge::end() const {
    return end_;
}

size_t Edge::common(const Sequence &other) const {
    size_t res = 0;
    while(res < seq.size() && res < other.size() && seq[res] == other[res]) {
        res += 1;
    }
    return res;
}

size_t Edge::size() const {
    return seq.size();
}

double Edge::getCoverage() const {
    return double(cov) / size();
}

size_t Edge::intCov() const {
    return cov;
}

void Edge::incCov(size_t val) const {
#pragma omp atomic
    cov += val;
}

bool Edge::operator==(const Edge &other) const {
    return this == &other;
}

bool Edge::operator!=(const Edge &other) const {
    return this != &other;
}

bool Edge::operator<(const Edge &other) const {
    if(this == &other)
        return false;
    if(start_ != other.start_)
        return *start_ < *other.start_;
    return this->seq < other.seq;
}

bool Edge::operator<=(const Edge &other) const {
    return *this == other || *this < other;
}

std::string Edge::str() const {
    std::stringstream ss;
    const dbg::Vertex &v = *start();
    ss << v.hash() << v.isCanonical() << "ACGT"[seq[0]];
    return ss.str();
}

Sequence Edge::kmerSeq(size_t pos) const {
    VERIFY(pos <= seq.size());
    size_t k = start_->seq.size();
    if (pos >= k)
        return seq.Subseq(pos - k, pos);
    else {
        return start()->seq.Subseq(pos) + seq.Subseq(0, pos);
    }
}

std::string Edge::getId() const {
    return start_->getId() + "ACGT"[seq[0]];
}

std::string Edge::getShortId() const {
    return start_->getShortId() + "ACGT"[seq[0]];
}

//std::ostream &operator<<(std::ostream &os, const Edge &edge) {
//    os << edge.getShortId();
//    return os;
//}

Vertex::Vertex(hashing::htype hash, Vertex *_rc) : hash_(hash), rc_(_rc), canonical(false) {
    omp_init_lock(&writelock);
}

bool Vertex::isCanonical() const {
    return canonical;
}

size_t Vertex::coverage() const {
    return coverage_;
}

bool Vertex::isCanonical(const Edge &edge) const {
    const Vertex &other = edge.end()->rc();
    if(hash() != other.hash())
        return hash() < other.hash();
    if (isCanonical() != other.isCanonical())
        return isCanonical();
    const Edge &rc_edge = edge.rc();
    return edge.seq <= rc_edge.seq;
}

std::string Vertex::edgeId(const Edge &edge) const {
    std::stringstream ss;
    if(isCanonical(edge)) {
        ss << hash() << isCanonical() << "ACGT"[edge.seq[0]];
        return ss.str();
    } else {
        return edge.end()->rc().edgeId(edge.rc());
    }
}

void Vertex::checkConsistency() const {
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

std::string Vertex::getId() const {
    std::stringstream ss;
    if(!isCanonical())
        ss << "-";
    ss << hash();
    return ss.str();
}

std::string Vertex::getShortId() const {
    std::stringstream ss;
    if(!isCanonical())
        ss << "-";
    ss << hash() % 10000000;
    return ss.str();
}

void Vertex::incCoverage() {
#pragma omp atomic update
    coverage_ += 1;
#pragma omp atomic update
    rc().coverage_ += 1;
}

void Vertex::setSequence(const Sequence &_seq) {
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

void Vertex::clearSequence() {
    if (!seq.empty()) {
        seq = Sequence();
        rc_->seq = Sequence();
    }
}

Edge &Vertex::addEdgeLockFree(const Edge &edge) {
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

void Vertex::addEdge(const Edge &e) {
    omp_set_lock(&writelock);
    addEdgeLockFree(e);
    omp_unset_lock(&writelock);
}

Edge &Vertex::getOutgoing(unsigned char c) const {
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

bool Vertex::hasOutgoing(unsigned char c) const {
    for (const Edge &edge : outgoing_) {
        if (edge.seq[0] == c) {
            return true;
        }
    }
    return false;
}

bool Vertex::operator<(const Vertex &other) const {
    return hash_ < other.hash_ || (hash_ == other.hash_ && canonical && !other.canonical);
}

void Vertex::clear() {
    outgoing_.clear();
    rc_->outgoing_.clear();
}

void Vertex::clearOutgoing() {
    outgoing_.clear();
}

Vertex::Vertex(hashing::htype hash) : hash_(hash), rc_(new Vertex(hash, this)), canonical(true) {
    omp_init_lock(&writelock);
}

Vertex::~Vertex() {
    if (rc_ != nullptr) {
        rc_->rc_ = nullptr;
        delete rc_;
    }
    rc_ = nullptr;
}

void Vertex::sortOutgoing() {
    std::sort(outgoing_.begin(), outgoing_.end());
}

bool Vertex::isJunction() const {
    return outDeg() != 1 || inDeg() != 1;
}

bool Vertex::operator==(const Vertex &other) const {
    return this == &other;
}

bool Vertex::operator!=(const Vertex &other) const {
    return this != &other;
}

void SparseDBG::checkSeqFilled(size_t threads, logging::Logger &logger) {
    logger.info() << "Checking vertex sequences" << std::endl;
    std::function<void(std::pair<const hashing::htype, Vertex> &)> task =
            [&logger](std::pair<const hashing::htype, Vertex> &pair) {
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

void SparseDBG::printEdge(std::ostream &os, Vertex &start, Edge &edge, bool output_coverage) {
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
