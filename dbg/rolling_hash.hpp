#pragma once
//
// Created by anton on 7/20/20.
//

template<typename T>
T pow(T base, size_t p) {
    if (p == 0)
        return 1;
    T tmp = pow(base, p / 2);
    if (p %2 == 1)
        return base * tmp * tmp;
    else
        return tmp * tmp;
}

const size_t w = 12000;

template<typename htype>
class RollingHash {
public:
    const size_t k;
    const size_t hbase;
    const htype kpow;

    RollingHash(size_t _k, htype _hbase) : k(_k), hbase(_hbase), kpow(pow(hbase, k - 1)){
    }

    htype hash(const Sequence &seq, size_t pos) const {
        htype hash = 0;
        for(size_t i = pos; i < pos + k; i++) {
            hash = hash * hbase + seq[i];
        }
        return hash;
    }

    htype shiftRight(const Sequence &seq, size_t pos, htype hash, char c) const {
        return (hash - kpow * dignucl(seq[pos])) * hbase + c;
    }

    htype next(const Sequence &seq, size_t pos, htype hash) const {
        return shiftRight(seq, pos, hash, seq.operator[](pos + k));
    }

    bool hasNext(const Sequence &seq, size_t pos) const {
        return pos + k < seq.size();
    }
};

template<typename htype>
class KWH {
public:
    const RollingHash<htype> & hasher;
    Sequence seq;
    size_t pos;
    htype hash;

    KWH(const RollingHash<htype> & _hasher, Sequence _seq, size_t _pos): hasher(_hasher), seq(_seq), pos(_pos), hash(_hasher.hash(_seq, _pos)) {
    }

    KWH(const RollingHash<htype> & _hasher, Sequence _seq, size_t _pos, htype _hash): hasher(_hasher), seq(_seq), pos(_pos), hash(_hash) {
    }

    KWH next() const {
        return {hasher, seq, pos + 1, hasher.next(seq, pos, hash)};
    }

    bool hasNext() const {
        return hasher.hasNext(seq, pos);
    }

    htype operator<<(char c) const {
        return hasher.shiftRight(c);
    }

    KWH &operator=(const KWH &other) {
        if(this == &other)
            return *this;
        seq = other.seq;
        pos = other.pos;
        hash = other.hash;
        return *this;
    }
};


template<typename htype>
class MinQueue {
    std::deque<KWH<htype>> q;
public:
    MinQueue() = default;

    void push(const KWH<htype> &kwh) {
        while(!q.empty() && q.back().hash >= kwh.hash) {
            q.pop_back();
        }
        q.push_back(kwh);
    }

    void pop(size_t pos) {
        if(!q.empty() && q.front().pos <= pos) {
            q.pop_front();
        }
    }

    bool empty() const {
        return q.empty();
    }

    KWH<htype> get() const {
        return q.front();
    }
};

template<typename htype>
class MinimizerCalculator {
private:
    const Sequence seq;
    const size_t w;
    KWH<htype> kwh;
    size_t pos;
    MinQueue<htype> queue;
public:
    MinimizerCalculator(const Sequence& _seq, const RollingHash<htype> &_hasher, size_t _w) :
            seq(_seq), w(_w), kwh(_hasher, seq, 0) {
        for(size_t i = 1; i + 1 < w; i++) {
            kwh = kwh.next();
            queue.push(kwh);
        }
    }

    htype next() {
        queue.pop(pos);
        pos += 1;
        kwh = kwh.next();
        queue.push(kwh);
        return queue.get().hash;
    }

    bool hasNext() const {
        return kwh.hasNext();
    }

    std::vector<htype> minimizers() {
        std::vector<htype> res;
        res.push_back(next());
        while(hasNext()) {
            htype val = next();
            if(val != res.back()) {
                res.push_back(val);
            }
        }
        return std::move(res);
    }
};

