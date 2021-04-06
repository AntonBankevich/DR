#pragma once

#include "ff.hpp"
#include "sparse_dbg.hpp"
#include "visualization.hpp"

class MappedNetwork : public Network {
public:
    std::unordered_map<int, dbg::Edge *> edge_mapping;
    std::unordered_map<dbg::Vertex *, int> vertex_mapping;
    MappedNetwork(const Component &component, size_t unique_len, size_t min_flow = 0) {
        dbg::SparseDBG &graph = component.graph;
        for(htype hash : component.v) {
            for(dbg::Vertex *v_it : {&graph.getVertex(hash), &graph.getVertex(hash).rc()}) {
                dbg::Vertex &v = *v_it;
                vertex_mapping[&v] = addVertex();
                std::cout << vertices.size() - 1 << " " << v.hash() << " " << v.isCanonical() << std::endl;
            }
        }
        for(htype hash : component.v) {
            for(dbg::Vertex *v_it : {&graph.getVertex(hash), &graph.getVertex(hash).rc()}) {
                dbg::Vertex &v = *v_it;
                for(dbg::Edge &edge : v.getOutgoing()) {
                    if(edge.size() < unique_len) {
                        int eid = addEdge(vertex_mapping[&v], vertex_mapping[edge.end()], min_flow, 10000);
                        edge_mapping[eid] = &edge;
                    } else {
                        addSink(vertex_mapping[&v], 1);
                        addSource(vertex_mapping[&v.rc()], 1);
                    }
                }
            }
        }
    }

};

class UniqueClassificator {
private:
    SparseDBG &dbg;
public:
    std::unordered_set<Edge *> unique_set;
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
            size_t min_flow = 1;
            MappedNetwork net(component, unique_len, min_flow);
            bool res = net.fillNetwork();
            if(res) {
                logger << "Found unique edges in component " << cnt << std::endl;
                std::unordered_map<int, size_t> multiplicities = net.findFixedMultiplicities();
                for (auto &rec : multiplicities) {
                    logger << "Edge " << net.edge_mapping[rec.first]->start()->hash()
                           << "ACGT"[net.edge_mapping[rec.first]->seq[0]]
                           << " has fixed multiplicity " << rec.second + min_flow << std::endl;
                    if (rec.second + min_flow == 1) {
                        unique_set.emplace(net.edge_mapping[rec.first]);
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
