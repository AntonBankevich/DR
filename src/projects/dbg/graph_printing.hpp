#pragma once

#include "component.hpp"
#include "sparse_dbg.hpp"
namespace dbg {
    inline void printFasta(std::ostream &out, const Component &component) {
        size_t cnt = 0;
        for (Edge &edge : component.edges()) {
            Sequence tmp = edge.start()->seq + edge.seq;
            Vertex &end = *edge.end();
            out << ">" << cnt << "_" << edge.start()->hash() << int(edge.start()->isCanonical()) <<
                "_" << end.hash() << int(end.isCanonical()) << "_" << edge.size() << "_" << edge.getCoverage()
                << std::endl;
            cnt++;
            out << tmp.str() << "\n";
        }
    }

    inline void printGFA(std::ostream &out, const Component &component, bool calculate_coverage) {
        out << "H\tVN:Z:1.0" << std::endl;
        size_t cnt = 0;
        for (Edge &edge : component.edges()) {
            if (edge.start()->isCanonical(edge)) {
                if (calculate_coverage)
                    out << "S\t" << edge.start()->edgeId(edge) << "\t" << edge.start()->seq << edge.seq
                        << "\tKC:i:" << edge.intCov() << std::endl;
                else
                    out << "S\t" << edge.start()->edgeId(edge) << "\t" << edge.start()->seq << edge.seq << std::endl;
            }
        }
        for (const auto &hash : component.v) {
            const Vertex &vertex = component.graph.getVertex(hash);
            for (const Edge &out_edge : vertex) {
                std::string outid = vertex.edgeId(out_edge);
                bool outsign = vertex.isCanonical(out_edge);
                for (const Edge &inc_edge : vertex.rc()) {
                    std::string incid = vertex.rc().edgeId(inc_edge);
                    bool incsign = !vertex.rc().isCanonical(inc_edge);
                    out << "L\t" << incid << "\t" << (incsign ? "+" : "-") << "\t" << outid << "\t"
                        << (outsign ? "+" : "-") << "\t" << component.graph.hasher().getK() << "M" << std::endl;
                }
            }
        }
    }
}