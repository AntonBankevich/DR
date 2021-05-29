#pragma once

#include <utility>


#include "ff.hpp"
#include "sparse_dbg.hpp"
#include "visualization.hpp"
#include "compact_path.hpp"

class MappedNetwork : public Network {
private:
    size_t min_flow;
public:
    std::unordered_map<int, dbg::Edge *> edge_mapping;
    std::unordered_map<dbg::Vertex *, int> vertex_mapping;

    MappedNetwork(const Component &component, const std::function<bool(const dbg::Edge &)> &unique, size_t min_flow = 0,
                  const std::function<size_t(const dbg::Edge &)> &max_flow = [](const dbg::Edge &){return 100000;});

    std::vector<dbg::Edge*> getUnique(logging::Logger &logger);
};

struct BoundRecord {
    size_t lowerBound;
    size_t upperBound;
    static size_t inf;
    BoundRecord() : lowerBound(0), upperBound(inf){
    }

    bool isUnique() const {
        return lowerBound == 1 && upperBound == 1;
    }
};

class MultiplicityBounds {
private:
    std::unordered_map<const Edge *, BoundRecord> multiplicity_bounds;
    size_t inf = 100000;
public:
    size_t upperBound(const Edge &edge) const {
        auto it = multiplicity_bounds.find(&edge);
        if(it == multiplicity_bounds.end())
            return inf;
        else return it->second.upperBound;
    }

    size_t lowerBound(const Edge &edge) const {
        auto it = multiplicity_bounds.find(&edge);
        if(it == multiplicity_bounds.end())
            return 0;
        else return it->second.lowerBound;
    }

    void setLowerBound(const Edge &edge, size_t val) {
        multiplicity_bounds[&edge].lowerBound = val;
    }

    void setUpperBound(const Edge &edge, size_t val) {
        multiplicity_bounds[&edge].upperBound = val;
    }

    void setBounds(const Edge &edge, size_t lower, size_t upper) {
        BoundRecord &bounds = multiplicity_bounds[&edge];
        bounds.lowerBound = lower;
        bounds.upperBound = upper;
    }

    bool isUnique(const Edge &edge) const {
        auto it = multiplicity_bounds.find(&edge);
        if(it == multiplicity_bounds.end())
            return false;
        return it->second.isUnique();
    }
};

class MultiplicityBoundsEstimator {
private:
    SparseDBG &dbg;
public:
    MultiplicityBounds bounds;

    void diploidUnique(size_t unique_length) {
    }

    void haploidUnique(size_t unique_length) {
    }

    std::vector<Component> uniqueSplit() {
        return {};
    }

    bool updateComponent(double rel_coverage) {
        return false;
    }

    bool updateComponentWithCoverage(double rel_coverage, double unique_coverage) {
        return false;
    }
};

class UniqueClassificator {
private:
    SparseDBG &dbg;
public:
    std::unordered_set<const Edge *> unique_set;
    const RecordStorage &reads_storage;
    explicit UniqueClassificator(SparseDBG &dbg, const RecordStorage &reads_storage) : dbg(dbg), reads_storage(reads_storage) {
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
            std::vector<const Edge *> new_unique = processComponent(logger, component, unique_len, out_file);
            unique_set.insert(new_unique.begin(), new_unique.end());
        }
    }

    std::vector<const dbg::Edge *> ProcessUsingCoverage(logging::Logger &logger, const Component &subcomponent,
                              const std::function<bool(const dbg::Edge &)> &is_unique, size_t min_flow) const {
        double max_cov = 0;
        double min_cov = 100000;
        for(htype hash : subcomponent.v) {
            for(Vertex *v_it : {&subcomponent.graph.getVertex(hash), &subcomponent.graph.getVertex(hash).rc()}) {
                Vertex &v = *v_it;
                for(Edge &edge : v) {
                    if(is_unique(edge)) {
                        const VertexRecord & record = reads_storage.getRecord(edge.end()->rc());
                        std::string s = edge.rc().seq.Subseq(0, 1).str();
                        size_t cnt = record.countStartsWith(Sequence(s + "A")) +
                                     record.countStartsWith(Sequence(s + "C")) +
                                     record.countStartsWith(Sequence(s + "G")) +
                                     record.countStartsWith(Sequence(s + "T"));
                        min_cov = std::min<double>(min_cov, cnt);
                        max_cov = std::max<double>(max_cov, cnt);
                    }
                }
            }
        }
        double threshold = std::max(min_cov * 1.4, max_cov * 1.1);
        std::function<size_t(const Edge &)> max_flow = [threshold](const Edge &edge) {
            return edge.getCoverage() < threshold ? 1 : 100000;
        };
        logger << "Attempting to use coverage for multiplicity estimation with coverage threshold " << threshold << std::endl;
        logger << "Component: ";
        for(htype hash : subcomponent.v) {
            logger << " " << hash % 100000;
        }
        logger << std::endl;
        MappedNetwork net2(subcomponent, is_unique, min_flow, max_flow);
        bool res2 = net2.fillNetwork();
        std::vector<const dbg::Edge *> extra_unique;
        if (res2) {
            logger << "Succeeded to use coverage for multiplicity estimation" << std::endl;
            for (Edge *edge : net2.getUnique(logger)) {
                extra_unique.emplace_back(edge);
            }
        } else {
            logger << "Failed to use coverage for multiplicity estimation" << std::endl;
        }
        return std::move(extra_unique);
    }

    std::vector<const dbg::Edge *> processComponent(logging::Logger &logger, const Component &component, size_t unique_len,
                                                    const std::experimental::filesystem::path &out_file) const {
        std::unordered_set<const dbg::Edge *> unique_in_component;
        size_t min_flow = 1;
        std::function<bool(const dbg::Edge &)> is_unique = [unique_len](const dbg::Edge &edge) {
            return edge.size() >= unique_len;
        };
        MappedNetwork net(component, is_unique, min_flow);
        bool res = net.fillNetwork();
        if(res) {
            logger << "Found unique edges in component" << std::endl;
            for(Edge * edge : net.getUnique(logger)) {
                unique_in_component.emplace(edge);
            }
        } else {
            logger << "Could not find unique edges in component" << std::endl;
            logger << "Relaxing flow conditions" << std::endl;
            min_flow = 0;
            MappedNetwork net1(component, is_unique, min_flow);
            bool res = net1.fillNetwork();
            if(res) {
                logger << "Found unique edges in component" << std::endl;
                for(Edge * edge : net1.getUnique(logger)) {
                    unique_in_component.emplace(edge);
                }
            } else {
                logger << "Could not find unique edges wth relaxed conditions in component" << std::endl;
            }
        }
        if(res) {
            std::function<bool(const dbg::Edge &)> is_unique = [&unique_in_component, unique_len](const dbg::Edge &edge){
                return edge.size() >= unique_len || unique_in_component.find(&edge) != unique_in_component.end();
            };
            std::vector<Component> subsplit = component.split(is_unique);
            for(Component &subcomponent : subsplit) {
                std::vector<const dbg::Edge *> extra_unique = ProcessUsingCoverage(logger, subcomponent, is_unique, min_flow);
                unique_in_component.insert(extra_unique.begin(), extra_unique.end());
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
        return {unique_in_component.begin(), unique_in_component.end()};
    }

    bool isUnique(const Edge &edge) const {
        return unique_set.find(&edge) != unique_set.end();
    }
};
