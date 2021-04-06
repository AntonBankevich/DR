#include "sparse_dbg.hpp"
using namespace dbg;

size_t Edge::updateTipSize() const {
    size_t new_val = 0;
    if(extraInfo == size_t(-1) && end_->inDeg() == 1) {
        for (const Edge & other : end_->getOutgoing()) {
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
    for(Edge &candidate : vend.getOutgoing()) {
        if(*candidate.end() == start_->rc() && candidate.size() == size() &&
           (size() <= seq.size() || candidate.seq.startsWith((!seq).Subseq(start_->seq.size())))) {
            return candidate;
        }
    }
    std::cout << start_->seq + seq << std::endl;
    std::cout << seq << std::endl;
    std::cout << vend.seq << std::endl;
    for(Edge &candidate : vend.getOutgoing()) {
        if(*candidate.end() == start_->rc() && candidate.size() == size() &&
           (size() <= seq.size() || candidate.seq.startsWith((!seq).Subseq(start_->seq.size())))) {
            std::cout << vend.seq + candidate.seq << std::endl;
        }
    }
    VERIFY(false);
    return vend.getOutgoing()[0];
}

Path Edge::walkForward() {
    Path res(*start_);
    res += *this;
    Vertex *next = end();
    while(next != nullptr && next != start_ && !next->isJunction()) {
        res += next->getOutgoing()[0];
        next = next->getOutgoing()[0].end();
    }
    return std::move(res);
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

bool Edge::operator<(const Edge &other) const {
    return this->seq < other.seq;
}


Vertex::Vertex(htype hash, Vertex *_rc) : hash_(hash), rc_(_rc), canonical(false) {
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