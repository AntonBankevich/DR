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

    inline void printFasta(const std::experimental::filesystem::path &outf, const Component &component) {
        std::ofstream out;
        out.open(outf);
        printFasta(out, component);
        out.close();
    }

    inline void printGFA(std::ostream &out, const Component &component, bool calculate_coverage) {
        out << "H\tVN:Z:1.0" << std::endl;
        size_t cnt = 0;
        for (Edge &edge : component.edges()) {
            if (edge.start()->isCanonical(edge)) {
                if (calculate_coverage)
                    out << "S\t" << edge.oldId() << "\t" << edge.start()->seq << edge.seq
                        << "\tKC:i:" << edge.intCov() << "\n";
                else
                    out << "S\t" << edge.oldId() << "\t" << edge.start()->seq << edge.seq << "\n";
            }
        }
        for (Vertex &vertex : component.verticesUnique()) {
            for (const Edge &out_edge : vertex) {
                std::string outid = out_edge.oldId();
                bool outsign = vertex.isCanonical(out_edge);
                for (const Edge &inc_edge : vertex.rc()) {
                    std::string incid = inc_edge.oldId();
                    bool incsign = !vertex.rc().isCanonical(inc_edge);
                    out << "L\t" << incid << "\t" << (incsign ? "+" : "-") << "\t" << outid << "\t"
                        << (outsign ? "+" : "-") << "\t" << component.graph().hasher().getK() << "M" << "\n";
                }
            }
        }
    }

    inline void printGFA(const std::experimental::filesystem::path &outf, const Component &component, bool calculate_coverage) {
        std::ofstream out;
        out.open(outf);
        printGFA(out, component, calculate_coverage);
        out.close();
    }
}