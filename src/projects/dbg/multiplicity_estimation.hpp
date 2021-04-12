#pragma once

#include <utility>


#include "ff.hpp"
#include "sparse_dbg.hpp"
#include "visualization.hpp"

class MappedNetwork : public Network {
private:
    size_t min_flow;
public:
    std::unordered_map<int, dbg::Edge *> edge_mapping;
    std::unordered_map<dbg::Vertex *, int> vertex_mapping;

    MappedNetwork(const Component &component, const std::function<bool(const dbg::Edge &)> &unique, size_t min_flow = 0,
                  const std::function<size_t(const dbg::Edge &)> &max_flow = [](const dbg::Edge &){return 100000;}) :
                min_flow(min_flow) {
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
                for(dbg::Edge &edge : v) {
                    if(unique(edge)) {
                        int eid = addEdge(vertex_mapping[&v], vertex_mapping[edge.end()], min_flow, max_flow(edge));
                        edge_mapping[eid] = &edge;
                    } else {
                        addSink(vertex_mapping[&v], 1);
                        addSource(vertex_mapping[&v.rc()], 1);
                    }
                }
            }
        }
    }

    std::vector<dbg::Edge*> getUnique(logging::Logger &logger) {
        std::vector<dbg::Edge*> res;
        std::unordered_map<int, size_t> multiplicities = findFixedMultiplicities();
        for (auto &rec : multiplicities) {
            logger << "Edge " << edge_mapping[rec.first]->start()->hash()
                   << "ACGT"[edge_mapping[rec.first]->seq[0]]
                   << " has fixed multiplicity " << rec.second + min_flow << std::endl;
            if (rec.second + min_flow == 1) {
                res.push_back(edge_mapping[rec.first]);
            }
        }
        return std::move(res);
    }
};

class UniqueClassificator {
private:
    SparseDBG &dbg;
public:
    std::unordered_set<Edge *> unique_set;
    explicit UniqueClassificator(SparseDBG &dbg) : dbg(dbg) {
    }

    void classify(logging::Logger &logger, size_t unique_len, const std::experimental::filesystem::path &dir) {
        Component graph(dbg);
        std::vector<Component> split = graph.split(unique_len);
        size_t cnt = 0;
        for(auto & it : dbg) {
            for(auto v_it : {&it.second, &it.second.rc()}) {
                for(Edge &edge : *v_it) {
                    if(edge.size() > unique_len) {
                        unique_set.emplace(&edge);
                    }
                }
            }
        }
        for(Component &component : split) {
            cnt += 1;
            std::experimental::filesystem::path out_file = dir / (std::to_string(cnt) + ".dot");
            std::vector<Edge *> new_unique = processComponent(logger, component, unique_len, out_file);
            unique_set.insert(new_unique.begin(), new_unique.end());
        }
    }

    std::vector<dbg::Edge *> processComponent(logging::Logger &logger, const Component &component, size_t unique_len,
                          const std::experimental::filesystem::path &out_file) const {
        std::unordered_set<Edge *> unique_in_component;
        size_t min_flow = 1;
        std::function<bool(const dbg::Edge &)> is_unique = [unique_len](const dbg::Edge &edge) {
            return edge.size() >= unique_len;
        };
        MappedNetwork net(component, is_unique, min_flow);
        bool res = net.fillNetwork();
        double max_cov = 0;
        for(htype hash : component.v) {
            for(dbg::Vertex *v_it : {&component.graph.getVertex(hash), &component.graph.getVertex(hash).rc()}) {
                dbg::Vertex &v = *v_it;
                for(dbg::Edge &edge : v) {
                    if(edge.size() >= unique_len) {
                        max_cov = std::max(max_cov, edge.getCoverage());
                    }
                }
            }
        }
        std::function<size_t(const dbg::Edge &)> max_flow = [max_cov](const dbg::Edge &edge) {
            return edge.getCoverage() < max_cov * 1.2 ? 1 : 100000;
        };
        if(res) {
            logger << "Found unique edges in component" << std::endl;
            for(Edge * edge : net.getUnique(logger)) {
                unique_in_component.emplace(edge);
            }
            logger << "Attempting to use coverage for multiplicity estimation" << std::endl;
            MappedNetwork net2(component, is_unique, min_flow, max_flow);
            bool res2 = net2.fillNetwork();
            if(res2) {
                logger << "Succeeded to use coverage for multiplicity estimation" << std::endl;
                for(Edge * edge : net2.getUnique(logger)) {
                    unique_in_component.emplace(edge);
                }
            } else {
                logger << "Failed to use coverage for multiplicity estimation" << std::endl;
            }
        } else {
            logger << "Could not find unique edges in component" << std::endl;
            logger << "Relaxing flow conditions" << std::endl;
            min_flow = 0;
            MappedNetwork net1(component, is_unique, min_flow);
            bool res1 = net1.fillNetwork();
            if(res1) {
                logger << "Found unique edges in component" << std::endl;
                for(Edge * edge : net1.getUnique(logger)) {
                    unique_in_component.emplace(edge);
                }
                logger << "Attempting to use coverage for multiplicity estimation" << std::endl;
                MappedNetwork net2(component, is_unique, min_flow, max_flow);
                bool res2 = net2.fillNetwork();
                if(res2) {
                    logger << "Succeeded to use coverage for multiplicity estimation" << std::endl;
                    for(Edge * edge : net2.getUnique(logger)) {
                        unique_in_component.emplace(edge);
                    }
                } else {
                    logger << "Failed to use coverage for multiplicity estimation" << std::endl;
                }
            } else {
                logger << "Could not find unique edges wth relaxed conditions in component" << std::endl;
            }
        }
        std::function<std::string(Edge &)> colorer = [&unique_in_component](Edge &edge) {
            return unique_in_component.find(&edge) == unique_in_component.end() ? "black" : "red";
        };
        logger << "Printing component to " << out_file << std::endl;
        const std::function<std::string(Edge &)> labeler = [](Edge &) {return "";};
        std::ofstream os;
        os.open(out_file);
        component.printDot(os, labeler, colorer);
        os.close();
        return {unique_set.begin(), unique_set.end()};
    }

    bool isUnique(Edge &edge) const {
        return unique_set.find(&edge) != unique_set.end();
    }
};
