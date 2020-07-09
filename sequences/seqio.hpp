#pragma once

#include <string>
#include "../common/string_utils.hpp"
#include "stream.hpp"
#include "contigs.hpp"

namespace io {

//    TODO: Deal with corrupted files, comments in read names
    class SeqReader {
    private:
        class Reader {
        public:
            virtual Contig read(std::istream &s) = 0;
        };

        class FastaReader: public Reader {
            Contig read(std::istream &s) override {
                string id, seq;
                std::getline(s, id);
                std::getline(s, seq);
                return Contig(Sequence(trim(seq)), trim(id.substr(1, id.size() - 1)));
            }
        };

        class FastqReader: public Reader {
            Contig read(std::istream &s) override {
                std::string id, seq;
                std::getline(s, id);
                std::getline(s, seq);
                Contig res(Sequence(trim(seq)), trim(id.substr(1, id.size() - 1)));
                std::getline(s, id);
                std::getline(s, seq);
                return res;
            }
        };

        std::istream * stream;
        Reader * reader;
    public:
        explicit SeqReader(const std::string & file_name) {
            if (endsWith(file_name, ".gz")) {
                stream = new gzstream::igzstream(file_name.c_str());
                if (endsWith(file_name, "fastq.gz") or endsWith(file_name, "fq.gz"))
                    reader = new FastqReader();
                else
                    reader = new FastaReader();
            } else {
                stream = new std::ifstream(file_name);
                if (endsWith(file_name, "fastq") or endsWith(file_name, "fq"))
                    reader = new FastqReader();
                else
                    reader = new FastaReader();
            }
        }

        Contig read() {
            return reader->read(*stream);
        }

        std::vector<Contig> readAll() {
            std::vector<Contig> res;
            while(!eof()) {
                res.emplace_back(reader->read(*stream));
            }
            return res;
        }

        bool eof() {
            return stream->eof();
        }

        ~SeqReader() {
            free(stream);
            free(reader);
        }
    };
}