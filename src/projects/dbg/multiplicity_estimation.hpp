#pragma once

#include "ff.hpp"
#include "sparse_dbg.hpp"
#include "visualization.hpp"

class UniqueClassificator {
private:
    SparseDBG &dbg;
    std::unordered_set<Edge *> unique_set;
public:
    UniqueClassificator(SparseDBG &dbg) : dbg(dbg) {
    }

    void classify(logging::Logger &logger, size_t unique_len, const std::experimental::filesystem::path &dir) {
        Component graph(dbg);
        std::vector<Component> split = graph.split(unique_len);
        size_t cnt = 0;
        for(auto & it : dbg) {
            for(auto v_it : {&it.second, &it.second.rc()}) {
                for(Edge &edge : v_it->getOutgoing()) {
                    if(edge.size() > unique_len) {
                        unique_set.emplace(&edge);
                    }
                }
            }
        }
        for(Component &component : split) {
            cnt += 1;
            std::unordered_map<int, Edge *> edge_mapping;
            std::unordered_map<Vertex *, int> vertex_mapping;
            Network net;
            for(htype hash : component.v) {
                for(Vertex *v_it : {&dbg.getVertex(hash), &dbg.getVertex(hash).rc()}) {
                    Vertex &v = *v_it;
                    vertex_mapping[&v] = net.addVertex();
                    std::cout << net.vertices.size() - 1 << " " << v.hash() << " " << v.isCanonical() << std::endl;
                }
            }
            for(htype hash : component.v) {
                for(Vertex *v_it : {&dbg.getVertex(hash), &dbg.getVertex(hash).rc()}) {
                    Vertex &v = *v_it;
                    for(Edge &edge : v.getOutgoing()) {
                        if(edge.size() < unique_len) {
                            int eid = net.addEdge(vertex_mapping[&v], vertex_mapping[edge.end()], 1, 10000);
                            edge_mapping[eid] = &edge;
                            std::cout << "New edge " << eid << " " << vertex_mapping[&v] << " " << vertex_mapping[edge.end()] << std::endl;
                        } else {
                            net.addSink(vertex_mapping[&v], 1);
                            net.addSource(vertex_mapping[&v.rc()], 1);
                            std::cout << "Sink " << vertex_mapping[&v] << std::endl;
                            std::cout << "Source " << vertex_mapping[&v.rc()] << std::endl;
                        }
                    }
                }
            }
            bool res = net.fillNetwork();
            if(res) {
                logger << "Found unique edges in component " << cnt << std::endl;
                std::unordered_map<int, size_t> multiplicities = net.findFixedMultiplicities();
                for (auto &rec : multiplicities) {
                    logger << "Edge " << edge_mapping[rec.first]->start()->hash()
                           << "ACGT"[edge_mapping[rec.first]->seq[0]]
                           << " has fixed multiplicity " << rec.second << std::endl;
                    if (rec.second == 0) {
                        unique_set.emplace(edge_mapping[rec.first]);
                    }
                }
            } else {
                logger << "Could not find unique edges in component " << cnt << std::endl;
            }
            std::function<std::string(Edge &)> colorer = [this](Edge &edge) {
                return unique_set.find(&edge) == unique_set.end() ? "black" : "red";
            };
            const std::function<std::string(Edge &)> labeler = [](Edge &) {return "";};
            std::ofstream os;
            os.open(dir / (std::to_string(cnt) + ".dot"));
            component.printDot(os, labeler, colorer);
            os.close();
        }
    }

    bool isUnique(Edge &edge) const {
        return unique_set.find(&edge) != unique_set.end();
    }
};
