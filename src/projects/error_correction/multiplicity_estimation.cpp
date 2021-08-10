#include "diploidy_analysis.hpp"
#include "multiplicity_estimation.hpp"

size_t BoundRecord::inf = 1000000000000ul;

MappedNetwork::MappedNetwork(const Component &component, const std::function<bool(const dbg::Edge &)> &unique,
                             double rel_coverage, double unique_coverage) {
    dbg::SparseDBG &graphr = component.graph();
    for(dbg::Vertex &v : component.vertices()) {
        vertex_mapping[&v] = addVertex();
    }
    for(dbg::Vertex &v : component.vertices()) {
        for(dbg::Edge &edge : v) {
            if (!unique(edge)) {
                size_t min_flow = (!edge.is_reliable || edge.getCoverage() < rel_coverage) ? 0 : 1;
                size_t max_flow = edge.getCoverage() < unique_coverage ? 1 : 1000000000;
                int eid = addEdge(vertex_mapping[&v], vertex_mapping[edge.end()], min_flow, max_flow);
                edge_mapping[eid] = &edge;
            } else {
                addSink(vertex_mapping[&v], 1);
                addSource(vertex_mapping[&v.rc()], 1);
            }
        }
    }
}

std::vector<dbg::Edge *> MappedNetwork::getUnique(logging::Logger &logger) {
    std::vector<dbg::Edge*> res;
    std::unordered_map<int, size_t> multiplicities = findFixedMultiplicities();
    for (auto &rec : multiplicities) {
//        logger << "Edge " << edge_mapping[rec.first]->start()->hash() << edge_mapping[rec.first]->start()->isCanonical()
//               << "ACGT"[edge_mapping[rec.first]->seq[0]]
//               << " has fixed multiplicity " << rec.second << std::endl;
        if(rec.second == 1)
            res.emplace_back(edge_mapping[rec.first]);
    }
    return std::move(res);
}

size_t MappedNetwork::addTipSinks() {
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

MultiplicityBoundsEstimator::MultiplicityBoundsEstimator(SparseDBG &dbg,
                                                         const AbstractUniquenessStorage &uniquenessStorage) : dbg(dbg){
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

bool MultiplicityBoundsEstimator::updateComponent(logging::Logger &logger, const Component &component,
                                                  const AbstractUniquenessStorage &uniquenessStorage,
                                                  double rel_coverage, double unique_coverage) {
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

void MultiplicityBoundsEstimator::update(logging::Logger &logger, double rel_coverage,
                                         const std::experimental::filesystem::path &dir) {
    ensure_dir_existance(dir);
    size_t cnt = 0;
    for(const Component &component : UniqueSplitter(bounds).splitGraph(dbg)) {
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

void UniqueClassificator::markPseudoHets() const {
    for(Edge &edge : dbg.edges()) {
        if(!isUnique(edge) || edge.end()->outDeg() != 2 || edge.end()->inDeg() != 1)
            continue;
        Vertex &start = *edge.end();
        Edge &correct = start[0].getCoverage() > start[1].getCoverage() ? start[0] : start[1];
        Edge &incorrect = start[0].getCoverage() <= start[1].getCoverage() ? start[0] : start[1];
        incorrect.is_reliable = false;
        incorrect.rc().is_reliable = false;
        GraphAlignment cor_ext = reads_storage.getRecord(start).
                getFullUniqueExtension(correct.seq.Subseq(0, 1), 1, 0).getAlignment();
        GraphAlignment incor_ext = reads_storage.getRecord(start).
                getFullUniqueExtension(incorrect.seq.Subseq(0, 1), 1, 0).getAlignment();
        for(Segment<Edge> &seg : incor_ext) {
            bool found= false;
            for(Segment<Edge> &seg1 : cor_ext) {
                if(seg.contig() == seg1.contig()) {
                    found = true;
                    break;
                }
            }
            if(found)
                break;
            seg.contig().is_reliable = false;
            seg.contig().rc().is_reliable = false;
        }
    }
}

void UniqueClassificator::classify(logging::Logger &logger, size_t unique_len,
                                   const std::experimental::filesystem::path &dir) {
    logger.info() << "Looking for unique edges" << std::endl;
    size_t cnt = 0;
    for(Edge &edge : dbg.edges()) {
        edge.is_reliable = true;
    }
    if(diploid) {
        SetUniquenessStorage duninque = BulgePathAnalyser(dbg, unique_len).uniqueEdges();
        for (Edge &edge : dbg.edges()) {
            if (duninque.isUnique(edge)) {
                addUnique(edge);
            }
        }
    } else {
        for (Edge &edge : dbg.edges()) {
            if(edge.size() > unique_len || (edge.start()->inDeg() == 0 && edge.size() > unique_len / 3)) {
                addUnique(edge);
            }
        }
    }
    markPseudoHets();
    std::vector<Component> split = UniqueSplitter(*this).split(Component(dbg));
    for(Component &component : split) {
        cnt += 1;
        std::experimental::filesystem::path out_file = dir / (std::to_string(cnt) + ".dot");
        printDot(out_file, component, reads_storage.labeler());
        logger.info() << "Component parameters: size=" << component.size() << " border=" << component.countBorderEdges() <<
                      " tips=" << component.countTips() <<
                      " subcomponents=" << component.realCC() << " acyclic=" << component.isAcyclic() <<std::endl;
        if(component.size() > 2 && component.countBorderEdges() == 2 &&component.countTips() == 0 &&
           component.realCC() == 2 && component.isAcyclic()) {
            processSimpleComponent(logger, component);
        }
        std::vector<const Edge *> new_unique = processComponent(logger, component);
        addUnique(new_unique.begin(), new_unique.end());
        logger << "Printing component to " << out_file << std::endl;
        const std::function<std::string(Edge &)> labeler = [](Edge &) {return "";};
        const std::function<std::string(Edge &)> colorer = [this](Edge &edge) {
            if(isUnique(edge)) {
                return "black";
            }
            if(!edge.is_reliable) {
                return "red";
            }
            return "blue";
        };
        printDot(out_file, component, reads_storage.labeler(), colorer);
    }
}

std::vector<const dbg::Edge *>
UniqueClassificator::ProcessUsingCoverage(logging::Logger &logger, const Component &subcomponent,
                                          const std::function<bool(const dbg::Edge &)> &is_unique,
                                          double rel_coverage) const {
    double max_cov = 0;
    double min_cov = 100000;
    for(Vertex &v : subcomponent.vertices()) {
        for(Edge &edge : v) {
            if(is_unique(edge)) {
                const VertexRecord & record = reads_storage.getRecord(edge.end()->rc());
                std::string s = edge.rc().seq.Subseq(0, 1).str();
                size_t cnt = record.countStartsWith(Sequence(s + "A")) +
                             record.countStartsWith(Sequence(s + "C")) +
                             record.countStartsWith(Sequence(s + "G")) +
                             record.countStartsWith(Sequence(s + "T"));
                min_cov = std::min<double>(min_cov, std::min<double>(cnt, edge.getCoverage()));
                max_cov = std::max<double>(max_cov, std::min<double>(cnt, edge.getCoverage()));
            }
        }
    }
    double threshold = std::max(min_cov * 1.4, max_cov * 1.2);
    double rel_threshold = std::min(min_cov * 0.9, max_cov * 0.7);
    logger << "Attempting to use coverage for multiplicity estimation with coverage threshold " << threshold << std::endl;
    logger << "Component: ";
    for(Vertex &vertex : subcomponent.verticesUnique()) {
        logger << " " << vertex.getShortId();
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

Edge &UniqueClassificator::getStart(const Component &component) const {
    for(Edge &edge : component.edges()) {
        if(!component.contains(*edge.end()))
            return edge.rc();
    }
    VERIFY(false);
}

void UniqueClassificator::processSimpleComponent(logging::Logger &logger, const Component &component) const {
    logger << "Collapsing acyclic component" << std::endl;
    typedef std::pair<size_t, Edge *> StoredValue;
    std::priority_queue<StoredValue> queue;
    queue.emplace(0, &getStart(component));
    std::unordered_map<Vertex *, Edge *> prev;
    Edge *end = nullptr;
    while(!queue.empty()) {
        auto tmp = queue.top();
        queue.pop();
        size_t score = tmp.first;
        Edge *last = tmp.second;
        Vertex *vert = last->end();
        if(prev.find(vert) != prev.end()) {
            continue;
        }
        prev[vert] = last;
        for(Edge &edge : *vert) {
            if(!component.contains(*edge.end())) {
                end = last;
                break;
            }
            queue.emplace(score + edge.intCov(), &edge);
        }
        if(end != nullptr)
            break;
    }
    if(end == nullptr) {
        logger << "Failed to collapse acyclic component" << std::endl;
        return;
    }
    for(Edge &edge : component.edgesInner()) {
        edge.is_reliable = false;
    }
    while(component.contains(*end->start())) {
        end->is_reliable = true;
        end = prev[end->start()];
    }
}

std::vector<const dbg::Edge *>
UniqueClassificator::processComponent(logging::Logger &logger, const Component &component) const {
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
        rel_coverage = 15;
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
        return (isUnique(edge) || unique_in_component.find(&edge) == unique_in_component.end()) ? "blue" : "black";
    };
    return {unique_in_component.begin(), unique_in_component.end()};
}
