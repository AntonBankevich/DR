#pragma once

#include <sequences/sequence.hpp>
#include <unordered_set>
#include <unordered_map>
#include <experimental/filesystem>
#include <fstream>
#include <common/string_utils.hpp>

namespace multigraph {
    class Edge;
    struct Vertex {
        Sequence seq;
        std::vector<Edge *> outgoing;
        Vertex *rc = nullptr;
        explicit Vertex(const Sequence &seq) : seq(seq) {
            outgoing.reserve(4);
        }

        size_t inDeg() const {
            return rc->outgoing.size();
        }

        size_t outDeg() const {
            return outgoing.size();
        }
    };

    struct Edge {
        Sequence seq;
        Vertex *start = nullptr;
        Vertex *end = nullptr;
        Edge *rc = nullptr;
        explicit Edge(const Sequence &seq) : seq(seq) {
        }
    };

    struct MultiGraph {
        std::vector<Vertex *> vertices;
        std::vector<Edge *> edges;
        MultiGraph() = default;
        MultiGraph(MultiGraph &&other) = default;
        MultiGraph(const MultiGraph &) = delete;

        ~MultiGraph() {
            for(Vertex * &v : vertices) {
                delete v;
                v = nullptr;
            }
            for(Edge * &e : edges) {
                delete e;
                e = nullptr;
            }
        }

        void checkConsistency() const {
            std::unordered_set<Edge const *> eset;
            std::unordered_set<Vertex const *> vset;
            for(Edge const * edge: edges) {
                eset.emplace(edge);
                VERIFY(edge->rc->start == edge->end->rc);
                VERIFY(edge->rc->rc == edge);
            }
            for(Edge const *edge : edges) {
                VERIFY(eset.find(edge->rc) != eset.end());
            }
            for(Vertex const * v: vertices) {
                vset.emplace(v);
                VERIFY(v->rc->rc == v);
                for(Edge *edge : v->outgoing) {
                    VERIFY(eset.find(edge) != eset.end());
                    VERIFY(edge->start == v);
                }
            }
            for(Vertex const *v : vertices) {
                VERIFY(vset.find(v->rc) != vset.end());
            }
        }

        Vertex &addVertex(const Sequence &seq) {
            vertices.emplace_back(new Vertex(seq));
            Vertex *res = vertices.back();
            Vertex *rc = res;
            if(seq != !seq) {
                vertices.emplace_back(new Vertex(!seq));
                rc = vertices.back();
            }
            res->rc = rc;
            rc->rc = res;
            return *res;
        }
        Edge &addEdge(Vertex &from, Vertex &to, const Sequence &seq) {
            edges.emplace_back(new Edge(seq));
            Edge *res = edges.back();
            res->start = &from;
            res->end = &to;
            res->start->outgoing.emplace_back(res);
            if(seq != !seq) {
                edges.emplace_back(new Edge(!seq));
                Edge *rc = edges.back();
                rc->start = to.rc;
                rc->end = from.rc;
                res->rc = rc;
                rc->rc = res;
                res->rc->start->outgoing.emplace_back(res->rc);
            } else {
                res->rc = res;
            }
            VERIFY(res->seq.startsWith(res->start->seq));
            VERIFY(res->rc->seq.startsWith(res->end->rc->seq));
            return *res;
        }

        std::vector<Edge *> uniquePathForward(Edge &edge) {
            std::vector<Edge *> res = {&edge};
            Vertex *cur = edge.end;
            while(cur != edge.start && cur->inDeg() == 1 && cur->outDeg() == 1) {
                res.emplace_back(cur->outgoing[0]);
                cur = res.back()->end;
            }
            return std::move(res);
        }

        std::vector<Edge *> uniquePath(Edge &edge) {
            std::vector<Edge *> path = uniquePathForward(*edge.rc);
            return uniquePathForward(*path.back()->rc);
        }

        MultiGraph Merge() {
            MultiGraph res;
            std::unordered_set<Edge *> used;
            std::unordered_map<Vertex *, Vertex *> old_to_new;
            for(Edge *edge : edges) {
                if(used.find(edge) != used.end())
                    continue;
                std::vector<Edge *> unique = uniquePath(*edge);
                for(Edge *e : unique) {
                    used.emplace(e);
                    used.emplace(e->rc);
                }
                Vertex *old_start = unique.front()->start;
                Vertex *old_end = unique.back()->end;
                Vertex *new_start = nullptr;
                Vertex *new_end = nullptr;
                if(old_to_new.find(old_start) == old_to_new.end()) {
                    old_to_new[old_start] = &res.addVertex(old_start->seq);
                    old_to_new[old_start->rc] = old_to_new[old_start]->rc;
                }
                new_start = old_to_new[old_start];
                if(old_to_new.find(old_end) == old_to_new.end()) {
                    old_to_new[old_end] = &res.addVertex(old_end->seq);
                    old_to_new[old_end->rc] = old_to_new[old_end]->rc;
                }
                new_end = old_to_new[old_end];
                Sequence new_seq;
                if(unique.size() == 1) {
                    new_seq = unique[0]->seq;
                } else {
                    SequenceBuilder sb;
                    sb.append(old_start->seq);
                    for(Edge *e : unique) {
                        sb.append(e->seq.Subseq(e->start->seq.size()));
                    }
                    new_seq = sb.BuildSequence();
                }
                res.addEdge(*new_start, *new_end, new_seq);
            }
            for(Vertex *vertex : vertices) {
                if(vertex->inDeg() == 0 && vertex->outDeg() == 0 && vertex->seq <= !vertex->seq) {
                    std::cout << "Isolated vertex " << vertex->seq.size() << " " << vertex->seq.Subseq(0, 20) << std::endl;
                    res.addVertex(vertex->seq);
                }
            }
            res.checkConsistency();
            return std::move(res);
        }

        void printCutEdges(const std::experimental::filesystem::path &f) {
            std::unordered_map<Vertex *, size_t> cut;
            for(Vertex *v : vertices) {
                if(v->seq <= !v->seq) {
                    if(v->outDeg() == 1) {
                        cut[v] = 0;
                    } else {
                        cut[v] = 1;
                    }
                    cut[v->rc] = 1 - cut[v];
                }
            }
            std::ofstream os;
            os.open(f);
            size_t cnt = 1;
            for(Edge *edge : edges) {
                if(edge->seq <= !edge->seq) {
                    size_t cut_left = edge->start->seq.size() * cut[edge->start];
                    size_t cut_right = edge->end->seq.size() * (1 - cut[edge->end]);
                    if(cut_left + cut_right + 1000 >= edge->seq.size()) {
                        continue;
                    }
                    os << ">" << cnt << "\n" << edge->seq.Subseq(cut_left, edge->seq.size() - cut_right) << "\n";
                }
            }
            os.close();
        }

        void printDot(const std::experimental::filesystem::path &f) {
            std::ofstream os;
            os.open(f);
            os << "digraph {\nnodesep = 0.5;\n";
            std::unordered_map<Vertex *, int> vmap;
            for(Vertex *vertex : vertices) {
                if(vertex->seq <= !vertex->seq) {
                    vmap[vertex] = vmap.size() + 1;
                    if(vertex->rc != vertex)
                        vmap[vertex->rc] = -vmap[vertex];
                }
            }
            for(Vertex *vertex : vertices) {
                os << vmap[vertex] << " [label=\"" << vertex->seq.size() << "\" style=filled fillcolor=\"white\"]\n";
            }
            size_t cnt = 1;
            std::unordered_map<Edge *, std::string> eids;
            for (Edge *edge : edges) {
                os << "\"" << vmap[edge->start] << "\" -> \"" << vmap[edge->end] <<
                   "\" [label=\"" << edge->seq.size() << "\" color = \"black\"]\n" ;
            }
            os << "}\n";
            os.close();
        }

        void printGFA(const std::experimental::filesystem::path &f) {
            std::ofstream os;
            os.open(f);
            os << "H\tVN:Z:1.0" << std::endl;
            size_t cnt = 1;
            std::unordered_map<Edge *, std::string> eids;
            for (Edge *edge : edges) {
                if (edge->seq <= !edge->seq) {
                    eids[edge] = itos(cnt);
                    eids[edge->rc] = itos(cnt);
                    cnt++;
                    os << "S\t" << eids[edge] << "\t" << edge->seq << "\n";
                }
            }
            for (Vertex *vertex : vertices) {
                if(!(vertex->seq <= !vertex->seq))
                    continue;
                for (Edge *out_edge : vertex->outgoing) {
                    std::string outid = eids[out_edge];
                    bool outsign = out_edge->seq <= !out_edge->seq;
                    for (Edge *inc_edge : vertex->rc->outgoing) {
                        std::string incid = eids[inc_edge];
                        bool incsign = !(out_edge->seq <= !out_edge->seq);
                        os << "L\t" << incid << "\t" << (incsign ? "+" : "-") << "\t" << outid << "\t"
                            << (outsign ? "+" : "-") << "\t" << vertex->seq.size() << "M" << "\n";
                    }
                }
            }
            os.close();
        }
    };
}