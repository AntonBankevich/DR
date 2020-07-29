//
// Created by anton on 7/27/20.
//
#pragma once


template<class htype>
class Edge;

template<class htype>
class Vertex {
private:
    std::vector<Edge> rightExtensions;
    std::vector<Edge> leftExtensions;
    omp_lock_t writelock{};
    Sequence seq;
    htype hash;
public:

    void clear() {
        rightExtensions.clear();
        leftExtensions.clear();
    }

    static void addSequence(const Sequence &s, std::vector<std::pair<Sequence, htype>> &extensions, htype next) {
        for(std::pair<Sequence, htype> & rightExtension : extensions) {
            if (s.size() <= rightExtension.first.size()) {
                if (s == rightExtension.first.Subseq(0, s.size())) {
                    return;
                }
            } else if (s.Subseq(0, rightExtension.first.size()) == rightExtension.first) {
                rightExtension = std::make_pair(Sequence(s.str()), next);
                return;
            }
        }
        extensions.emplace_back(Sequence(s.str()), next);
    }
    lock() {
        omp_set_lock(&writelock);
    }

    unlock() {
        omp_unset_lock(&writelock);
    }

public:
    friend class AssemblyGraph<htype>;
    ExtensionsRecord() {
        omp_init_lock(&writelock);
    }

    void setSequence(const Sequence &_seq) {
        if(seq.empty()) {
            seq = Sequence(_seq.str());
        }
    }

    void addRightExtension(const Sequence &s) {
        omp_set_lock(&writelock);
        addSequence(s, rightExtensions);
        omp_unset_lock(&writelock);
    }

    void addLeftExtension(const Sequence &s) {
        omp_set_lock(&writelock);
        addSequence(s, leftExtensions);
        omp_unset_lock(&writelock);
    }

    void addExtension(const Sequence &left, const Sequence &right, htype prev, htype next) {
        omp_set_lock(&writelock);
        addSequence(left, rightExtensions, next);
        addSequence(right, leftExtensions, prev);
        omp_unset_lock(&writelock);
    }

    std::vector<std::pair<Sequence, htype>> &getRightExtensions() {
        return rightExtensions;
    }
    std::vector<std::pair<Sequence, htype>> &getLeftExtensions() {
        return leftExtensions;
    }

    bool isJunction() const {
        return leftExtensions.size() != 1 || rightExtensions.size() != 1;
    }
};

template<class htype>
class Edge {
private:
    Sequence seq;
    htype otherVertex;
public:
    Edge(const Sequence &_s, type _otherVertex) : seq(_s), otherVertex(_otherVertex){
    }
};


template<typename htype>
class AssemblyGraph {
private:
public:
//    TODO: replace with perfect hash map? It is parallel, maybe faster and compact.
    std::unordered_map<htype, ExtensionsRecord<htype>> v;
    const RollingHash<htype> hasher;
    size_t w;
    SparseDBG(std::vector<htype> &&_hashs, RollingHash<htype> _hasher, size_t _w) : hasher(_hasher), w(_w) {
        std::vector<htype> hashs(_hashs);
        for(htype hash: hashs) {
            v[hash] = ExtensionsRecord<htype>();
        }
    }

    SparseDBG(SparseDBG<htype> &&other): hasher(other.hasher) {
        std::swap(v, other.v);
        w = other.w;
    }

    void addVertex(htype h) {
        v[h] = ExtensionsRecord<htype>();
    }

    std::vector<KWH<htype>> extractMinimizers(const Sequence &seq) const {
        std::vector<KWH<htype>> res;
        KWH<htype> kwh(hasher, seq, 0);
        while(true) {
            if(v.find(kwh.hash) != v.end()) {
                res.emplace_back(kwh);
            }
            if (!kwh.hasNext())
                break;
            kwh = kwh.next();
        }
        return std::move(res);
    }

    void processRead(const Sequence & _seq) {
        Sequence seq(_seq.str());
        std::vector<KWH<htype>> kmers = extractMinimizers(seq);
        std::vector<std::pair<size_t, htype>> pos;
        pos.push_back(std::make_pair(0, htype(-1)));
        for(size_t i = 0; i < kmers.size(); i++) {
            pos.push_back(std::make_pair(kmers[i].pos, kmers[i].hash));
        }
        pos.push_back(std::make_pair(seq.size() - hasher.k, htype(-1)));
        for(size_t i = 1; i + 1 < pos.size(); i++) {
            ExtensionsRecord<htype> &rec = v[pos[i].second];
//            TODO: if too memory heavy save only some of the labels
            rec.setSequence(seq.Subseq(pos[i].first, pos[i].first + hasher.k));
            rec.addExtension(!seq.Subseq(pos[i - 1].first, pos[i].first),
                             seq.Subseq(pos[i].first + hasher.k, pos[i + 1].first + hasher.k),
                             pos[i - 1].second, pos[i + 1].second);
        }
    }

    Sequence walkForward(htype start_hash, ExtensionsRecord<htype> &rec, size_t t) {
        std::vector<Sequence> res;
        if(rec.getLeftExtensions().size() > 0) {
            res.push_back(!(rec.getLeftExtensions()[0].first));
        }
        res.push_back(rec.seq);
        res.push_back(rec.getRightExtensions()[t].first);
        htype hash = rec.getRightExtensions()[t].second;
        while(hash != htype(-1) && hash != start_hash) {
            ExtensionsRecord<htype> &next = v.find(hash)->second;
            if (next.isJunction()) {
                if(next.getRightExtensions().size() > 0) {
                    res.push_back(next.getRightExtensions()[0].first);
                }
                break;
            }
            res.push_back(next.getRightExtensions()[0].first);
            hash = next.getRightExtensions()[0].second;
        }
        return Sequence::Concat(res);
    }

    void printStats() {
        std::vector<size_t> arr(10);
        size_t n11 = 0;
        size_t n01 = 0;
        for(auto &val: v) {
            ExtensionsRecord<htype> &tmp = val.second;
            arr[std::min(arr.size() - 1, tmp.getLeftExtensions().size())] += 1;
            arr[std::min(arr.size() - 1, tmp.getRightExtensions().size())] += 1;
            if (tmp.getLeftExtensions().size() == 1 && tmp.getRightExtensions().size() == 1) {
                n11 += 1;
            }
            if (tmp.getLeftExtensions().size() + tmp.getRightExtensions().size() == 1) {
                n01 += 1;
            }
        }
        size_t total = std::accumulate(arr.begin(), arr.end(), size_t(0));
        std::cout << "Sparse graph statistics:" << std::endl << "Total edges: " << total << std::endl;
        std::cout << "Number of end vertices: " << n01 << std::endl;
        std::cout << "Number of unbranching vertices: " << n11 << std::endl << "Distribution of degrees:" << std::endl;
        for(size_t i = 0; i < arr.size(); i++) {
            std::cout << i << " " << arr[i] << std::endl;
        }
    }

};
