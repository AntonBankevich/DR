#pragma once

#include <string>
#include "../common/string_utils.hpp"
#include "../common/stream.hpp"
#include "contigs.hpp"

namespace io {

    class SeqReader {
    private:
        class Reader {
        public:
            virtual Contig read(std::istream &s) = 0;
        };

        class FastaReader: public Reader {
            Contig read(std::istream &s) {
                string id, seq;
                std::getline(s, id);
                std::getline(s, seq);
                return Contig(seq, id);
            }
        };

        std::basic_istream<char> * stream;
    public:
        SeqReader(const std::string & file_name) {
            if (endsWith(file_name, ".gz")) {
                stream = new gzstream::igzstream(file_name.c_str());
            } else {
                stream = new std::ifstream(file_name);
            }
        }

        ~SeqReader() {
            free(stream);
        }
    };
}