#include "multiplicity_estimation.hpp"

size_t BoundRecord::inf = 1000000000000ul;

MappedNetwork::MappedNetwork(const Component &component, const std::function<bool(const dbg::Edge &)> &unique,
                             size_t min_flow, const std::function<size_t(const dbg::Edge &)> &max_flow) :
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
                if(!unique(edge)) {
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

std::vector<dbg::Edge *> MappedNetwork::getUnique(logging::Logger &logger) {
    std::vector<dbg::Edge*> res;
    std::unordered_map<int, size_t> multiplicities = findFixedMultiplicities();
    for (auto &rec : multiplicities) {
        logger << "Edge " << edge_mapping[rec.first]->start()->hash() << edge_mapping[rec.first]->start()->isCanonical()
               << "ACGT"[edge_mapping[rec.first]->seq[0]]
               << " has fixed multiplicity " << rec.second + min_flow << std::endl;
        if (rec.second + min_flow == 1) {
            res.push_back(edge_mapping[rec.first]);
        }
    }
    return std::move(res);
}
