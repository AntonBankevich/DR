#pragma once

#include <utility>


#include "ff.hpp"
#include "sparse_dbg.hpp"
#include "visualization.hpp"
#include "compact_path.hpp"

class MappedNetwork : public Network {
public:
    std::unordered_map<int, dbg::Edge *> edge_mapping;
    std::unordered_map<dbg::Vertex *, int> vertex_mapping;

    MappedNetwork(const Component &component, const std::function<bool(const dbg::Edge &)> &unique, double rel_coverage = 1000, double unique_coverage = 0);

    size_t addTipSinks() {
        size_t res = 0;
        for(Edge &edge : edges) {
            if(edge.min_flow == 0)
                continue;
            Vertex end = vertices[edge.end];
            if(end.out.empty()) {
                addSink(end.id, 1);
                res += 1;
            }
            Vertex start = vertices[edge.start];
            if(start.out.empty()) {
                addSource(start.id, 1);
                res += 1;
            }
        }
        return res;
    }

    std::vector<dbg::Edge*> getUnique(logging::Logger &logger);
};

class AbstractUniquenessStorage {
private:
    bool checkEdge(const dbg::Edge &edge) const {
        if(edge.end()->outDeg() + 1 != edge.end()->inDeg())
            return false;
        for(const Edge &e : *edge.end()) {
            if(!isUnique(e))
                return false;
        }
        for(const Edge &e : edge.end()->rc()) {
            if(e != edge.rc() && !isUnique(e))
                return false;
        }
        return true;
    }
public:
    virtual bool isUnique(const dbg::Edge &) const = 0;
    virtual ~AbstractUniquenessStorage() = default;

    bool isError(const dbg::Edge &edge) const {
        if(isUnique(edge))
            return false;
        return checkEdge(edge) || checkEdge(edge.rc());
    }

    std::function<std::string(const Edge &)> colorer(const std::string &unique_color = "black",
                                                     const std::string &repeat_color = "blue") const {
        return [this, unique_color, repeat_color](const Edge &edge) -> std::string {
            if(isUnique(edge))
                return unique_color;
            else
                return repeat_color;
        };
    }
};

class UniqueSplitter : public ConditionSplitter {
public:
    explicit UniqueSplitter(const AbstractUniquenessStorage &storage) :
            ConditionSplitter([&storage](const Edge& edge){return storage.isUnique(edge);}){
    }
};


class SetUniquenessStorage : public AbstractUniquenessStorage{
private:
    std::unordered_set<const dbg::Edge *> unique;
public:
    SetUniquenessStorage() {
    }

    template<class I>
    SetUniquenessStorage(I begin, I end) {
        addUnique(begin, end);
    }

    bool isUnique(const dbg::Edge &edge) const override {
        return unique.find(&edge) != unique.end();
    }

    void addUnique(const dbg::Edge &edge) {
        unique.emplace(&edge);
        unique.emplace(&edge.rc());
    }

    template<class I>
    void addUnique(I begin, I end) {
        while(begin != end) {
            const dbg::Edge &edge = **begin;
            unique.emplace(&edge);
            unique.emplace(&edge.rc());
            ++begin;
        }
    }
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

    size_t updateLowerBound(size_t val) {
        lowerBound = std::max(lowerBound, val);
        return lowerBound;
    }

    size_t updateUpperBound(size_t val) {
        upperBound = std::min(upperBound, val);
        return upperBound;
    }
};

class MultiplicityBounds : public AbstractUniquenessStorage {
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

    void updateLowerBound(const Edge &edge, size_t val) {
        BoundRecord &bounds = multiplicity_bounds[&edge];
        bounds.updateLowerBound(val);
    }

    void updateUpperBound(const Edge &edge, size_t val) {
        BoundRecord &bounds = multiplicity_bounds[&edge];
        bounds.updateUpperBound(val);
    }

    void updateBounds(const Edge &edge, size_t lower, size_t upper) {
        BoundRecord &bounds = multiplicity_bounds[&edge];
        bounds.updateLowerBound(lower);
        bounds.updateUpperBound(upper);
    }

    bool isUnique(const Edge &edge) const override {
        auto it = multiplicity_bounds.find(&edge);
        if(it == multiplicity_bounds.end())
            return false;
        return it->second.isUnique();
    }

    std::function<std::string(const Edge &)> labeler() const {
        return [this](const Edge &edge) -> std::string {
            auto it = multiplicity_bounds.find(&edge);
            if(it == multiplicity_bounds.end()) {
                return "";
            }
            std::stringstream ss;
            ss << "[" << it->second.lowerBound << "-";
            size_t upper = it->second.upperBound;
            if(upper < 100)
                ss << upper;
            else
                ss << "inf";
            ss << "]";
            return ss.str();
        };
    }
};

class MultiplicityBoundsEstimator {
private:
    SparseDBG &dbg;
    MultiplicityBounds bounds;
public:
    MultiplicityBoundsEstimator(SparseDBG &dbg, const AbstractUniquenessStorage &uniquenessStorage) : dbg(dbg){
        for(auto & it : dbg) {
            for(auto v_it : {&it.second, &it.second.rc()}) {
                for(Edge &edge : *v_it) {
                    if (uniquenessStorage.isUnique(edge)) {
                        bounds.updateBounds(edge, 1, 1);
                    }
                }
            }
        }
    }

    bool updateComponent(logging::Logger &logger, const Component &component, const AbstractUniquenessStorage &uniquenessStorage,
                                double rel_coverage, double unique_coverage = 0) {
        std::unordered_set<const dbg::Edge *> unique_in_component;
        std::function<bool(const dbg::Edge &)> is_unique =
                [&uniquenessStorage, unique_coverage, rel_coverage](const dbg::Edge &edge) {
            return uniquenessStorage.isUnique(edge) || (edge.getCoverage() >= rel_coverage && edge.getCoverage() < unique_coverage);
        };
        MappedNetwork net(component, is_unique, rel_coverage);
        bool res = net.fillNetwork();
        if(!res) {
            logger << "Initial flow search failed. Adding tip sinks." << std::endl;
            size_t tips = net.addTipSinks();
            if(tips > 0)
                res = net.fillNetwork();
        }
        if(res) {
            logger << "Found multiplicity bounds in component" << std::endl;
            for(auto rec : net.findBounds()) {
                this->bounds.updateBounds(*net.edge_mapping[rec.first], rec.second.first, rec.second.second);
            }
            return true;
        }
        logger << "Flow search failed. Multiplicity bounds were not updated." << std::endl;
        return false;
    }

    void update(logging::Logger &logger, double rel_coverage, const std::experimental::filesystem::path &dir) {
        ensure_dir_existance(dir);
        size_t cnt = 0;
        for(const Component &component : UniqueSplitter(bounds).split(dbg)) {
            if(component.size() <= 2)
                continue;
            cnt += 1;
            std::ofstream os1;
            os1.open(dir / (std::to_string(cnt) + "_before.dot"));
            printDot(os1, component, bounds.labeler(), bounds.colorer());
            os1.close();
            updateComponent(logger, component, bounds, rel_coverage);
            std::experimental::filesystem::path out_file = dir / (std::to_string(cnt) + ".dot");
            logger << "Printing component to " << out_file << std::endl;
            std::ofstream os;
            os.open(out_file);
            printDot(os, component, bounds.labeler(), bounds.colorer());
            os.close();
        }
    }
};

class UniqueClassificator : public SetUniquenessStorage{
private:
    SparseDBG &dbg;
public:
    const RecordStorage &reads_storage;
    void classify(logging::Logger &logger, size_t unique_len, const std::experimental::filesystem::path &dir) {
        size_t cnt = 0;
        for(Edge &edge : dbg.edges()) {
            edge.is_reliable = true;
        }
        for(Edge &edge : dbg.edges()) {
            if(edge.size() > unique_len || (edge.start()->inDeg() == 0 && edge.size() > unique_len / 3)) {
                addUnique(edge);
                if(edge.end()->outDeg() > 1 && edge.end()->inDeg() == 1) {
                    Edge *best = nullptr;
                    for(Edge &out : *edge.end()) {
                        if(best == nullptr) {
                            best = &out;
                            continue;
                        } else {
                            if (out.getCoverage() > best->getCoverage()) {
                                best->is_reliable = false;
                                best = &out;
                            } else {
                                out.is_reliable = false;
                            }
                        }
                    }
                }
            }
        }
        std::vector<Component> split = UniqueSplitter(*this).split(Component(dbg));
        for(Component &component : split) {
            cnt += 1;
            std::experimental::filesystem::path out_file = dir / (std::to_string(cnt) + ".dot");
            std::vector<const Edge *> new_unique = processComponent(logger, component, out_file);
            addUnique(new_unique.begin(), new_unique.end());
        }
    }

    explicit UniqueClassificator(SparseDBG &dbg, const RecordStorage &reads_storage) : dbg(dbg), reads_storage(reads_storage) {
    }

    std::vector<const dbg::Edge *> ProcessUsingCoverage(logging::Logger &logger, const Component &subcomponent,
                              const std::function<bool(const dbg::Edge &)> &is_unique, double rel_coverage) const {
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
        double threshold = std::max(min_cov * 1.4, max_cov * 1.2);
        double rel_threshold = std::min(min_cov * 0.9, max_cov * 0.7);
        logger << "Attempting to use coverage for multiplicity estimation with coverage threshold " << threshold << std::endl;
        logger << "Component: ";
        for(htype hash : subcomponent.v) {
            logger << " " << hash % 100000;
        }
        logger << std::endl;
        MappedNetwork net2(subcomponent, is_unique, rel_coverage, threshold);
        bool res2 = net2.fillNetwork();
        std::vector<const dbg::Edge *> extra_unique;
        if (res2) {
            logger << "Succeeded to use coverage for multiplicity estimation" << std::endl;
            for (Edge *edge : net2.getUnique(logger)) {
                extra_unique.emplace_back(edge);
            }
        } else {
            logger << "Failed to use coverage for multiplicity estimation" << std::endl;
            if(rel_coverage == 0) {
                logger << "Adjusted reliable edge threshold from " << rel_coverage << " to " << rel_threshold << std::endl;
                MappedNetwork net3(subcomponent, is_unique, rel_threshold, threshold);
                bool res3 = net3.fillNetwork();
                if (res3) {
                    logger << "Succeeded to use coverage for multiplicity estimation" << std::endl;
                    for (Edge *edge : net3.getUnique(logger)) {
                        extra_unique.emplace_back(edge);
                    }
                } else {
                    logger << "Failed to use coverage for multiplicity estimation" << std::endl;
                }
            }
        }
        return std::move(extra_unique);
    }

    std::vector<const dbg::Edge *> processComponent(logging::Logger &logger, const Component &component,
                                                    const std::experimental::filesystem::path &out_file) const {
        std::unordered_set<const dbg::Edge *> unique_in_component;
        double rel_coverage = 0;
        std::function<bool(const dbg::Edge &)> is_unique = [this](const dbg::Edge &edge) {
            return isUnique(edge);
        };
        MappedNetwork net(component, is_unique, rel_coverage);
        bool res = net.fillNetwork();
        if(res) {
            logger << "Found unique edges in component" << std::endl;
            for(Edge * edge : net.getUnique(logger)) {
                unique_in_component.emplace(edge);
            }
        } else {
            logger << "Could not find unique edges in component" << std::endl;
            logger << "Relaxing flow conditions" << std::endl;
            rel_coverage = 20;
            MappedNetwork net1(component, is_unique, rel_coverage);
            res = net1.fillNetwork();
            if(res) {
                logger << "Found unique edges in component" << std::endl;
                for(Edge * edge : net1.getUnique(logger)) {
                    unique_in_component.emplace(edge);
                }
            } else {
                logger << "Could not find unique edges with relaxed conditions in component" << std::endl;
            }
        }
        if(res) {
            std::function<bool(const dbg::Edge &)> is_unique1 = [&unique_in_component, this](const dbg::Edge &edge){
                return isUnique(edge) || unique_in_component.find(&edge) != unique_in_component.end();
            };
            std::vector<Component> subsplit = ConditionSplitter(is_unique1).split(component);
            for(Component &subcomponent : subsplit) {
                std::vector<const dbg::Edge *> extra_unique = ProcessUsingCoverage(logger, subcomponent, is_unique1, rel_coverage);
                unique_in_component.insert(extra_unique.begin(), extra_unique.end());
            }
        }
        std::function<std::string(Edge &)> colorer = [this, &unique_in_component](Edge &edge) {
            return (isUnique(edge) || unique_in_component.find(&edge) == unique_in_component.end()) ? "black" : "red";
        };
        logger << "Printing component to " << out_file << std::endl;
        const std::function<std::string(Edge &)> labeler = [](Edge &) {return "";};
        std::ofstream os;
        os.open(out_file);
        printDot(os, component, labeler, colorer);
        os.close();
        return {unique_in_component.begin(), unique_in_component.end()};
    }
};
