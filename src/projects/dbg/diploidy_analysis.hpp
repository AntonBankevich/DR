#pragma once

#include "sparse_dbg.hpp"
#include "multiplicity_estimation.hpp"

class BulgePath {
private:
    typedef typename std::vector<std::pair<dbg::Edge *, dbg::Edge *>> storage_type;
    typedef typename storage_type::const_iterator iterator_type;
    dbg::Vertex *start_;
    std::vector<std::pair<dbg::Edge *, dbg::Edge *>> path;
public:
    BulgePath(dbg::Vertex &_start) : start_(&_start) {
    }

    BulgePath(dbg::Edge &edge) : start_(edge.start()) {
        path.emplace_back(&edge, &edge);
    }

    BulgePath(std::vector<std::pair<dbg::Edge *, dbg::Edge *>> &&path_) : path(path_), start_(path[0].first->start()) {
    }

    dbg::Vertex &finish() const {
        if(path.empty())
            return *start_;
        return *path.back().first->end();
    }

    dbg::Vertex &start() const {
        return *start_;
    }

    void extend() {
        dbg::Vertex &last = finish();
        size_t deg = last.outDeg();
        VERIFY(deg == 1 || deg == 2);
        path.emplace_back(&last[0], &last[deg - 1]);
    }

    BulgePath RC() {
        if(path.empty()) {
            return BulgePath(start_->rc());
        }
        std::vector<std::pair<dbg::Edge *, dbg::Edge *>> rc;
        for(size_t i = 0; i < path.size(); i++) {
            rc.emplace_back(&path[path.size() - 1 - i].first->rc(), &path[path.size() - 1 - i].second->rc());
        }
        return BulgePath(std::move(rc));
    }

    BulgePath operator+(const BulgePath &other) const {
        VERIFY(finish() == other.start())
        std::vector<std::pair<dbg::Edge *, dbg::Edge *>> sum(path);
        sum.insert(sum.end(), other.path.begin(), other.path.end());
        return {std::move(sum)};
    }

    const std::pair<dbg::Edge *, dbg::Edge *> &operator[](size_t ind) const {
        return path[ind];
    }

    iterator_type begin() const {
        return path.begin();
    }

    iterator_type end() const {
        return path.end();
    }

    dbg::Vertex &vertexAt(size_t ind) {
        if(ind == 0)
            return *start_;
        return *path[ind - 1].first->end();
    }

    size_t size() const {
        return path.size();
    }

    size_t length() const {
        size_t res = 0;
        for(auto & p : path) {
            res += std::max(p.first->size(), p.second->size());
        }
        return res;
    }

    std::string str() const {
        std::stringstream ss;
        ss << start().label();
        for(const auto &p : path) {
            if(p.first == p.second) {
                ss << "-" << p.first->size() << "ACGT"[p.first->seq[0]] << "-" << p.first->end()->label();
            } else {
                ss << "-(" << p.first->size() << "ACGT"[p.first->seq[0]] << "," <<
                        p.second->size() << "ACGT"[p.second->seq[0]] << ")-" << p.first->end()->label();
            }
        }
        return ss.str();
    }
};


class BulgePathAnalyser {
private:

    static bool checkVertex(const dbg::Vertex &v) {
        return v.outDeg() == 1 ||
               (v.outDeg() == 2 && v[0].end() == v[1].end() && v[0].size() < v[1].size() * 1.3 && v[1].size() < v[0].size() * 1.3);
    }

    static bool isInner(const dbg::Vertex &v) {
        return checkVertex(v) && checkVertex(v.rc());
    }

    static BulgePath forwardPath(dbg::Vertex &start) {
        BulgePath res(start);
        dbg::Vertex * cur = &start;
        while(checkVertex(*cur)) {
            res.extend();
            cur = &res.finish();
            if(cur == &start)
                return std::move(res);
        }
        return std::move(res);
    }

    dbg::SparseDBG &dbg;
public:
    std::vector<BulgePath> paths;

    BulgePathAnalyser(dbg::SparseDBG &dbg, size_t min_len = 100000) : dbg(dbg) {
        std::unordered_set<dbg::Vertex *> visited;
        for(auto &it : dbg) {
            if(visited.find(&it.second) != visited.end())
                continue;
            if(checkVertex(it.second)) {
                BulgePath new_path = forwardPath(it.second);
                VERIFY(new_path.size() > 0);
                if(new_path.start() != new_path.finish()) {
                    BulgePath p2 = forwardPath(it.second.rc());
                    BulgePath p3 = p2.RC();
                    new_path = p3 + new_path;
                }
                if(new_path.length() < min_len)
                    continue;
                for(size_t i = 1; i + 1 <= new_path.size(); i++) {
                    visited.emplace(&new_path.vertexAt(i));
                    visited.emplace(&new_path.vertexAt(i).rc());
                }
                if(new_path.start() == new_path.finish()) {
                    visited.emplace(&new_path.start());
                    visited.emplace(&new_path.start().rc());
                }
                paths.emplace_back(new_path.RC());
                paths.emplace_back(std::move(new_path));
            }
        }
        for(auto &it : dbg) {
            if(visited.find(&it.second) != visited.end())
                continue;
            for(dbg::Vertex * vit : {&it.second, &it.second.rc()}) {
                for(dbg::Edge & edge : *vit) {
                    if(edge.size() > min_len) {
                        paths.emplace_back(edge);
                    }
                }
            }
        }
    }

    SetUniquenessStorage uniqueEdges() {
        std::vector<dbg::Edge *> res;
        for(BulgePath &bp : paths) {
            for (auto & p : bp) {
                if(p.first != p.second) {
                    res.emplace_back(p.first);
                    res.emplace_back(p.second);
                }
            }
            if(bp.size() == 1)
                res.emplace_back(bp[0].first);
        }
        return {res.begin(), res.end()};
    }
};