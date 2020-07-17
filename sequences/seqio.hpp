#pragma once

#include <string>
#include "../common/string_utils.hpp"
#include "stream.hpp"
#include "contigs.hpp"

namespace io {

//    TODO: Deal with corrupted files, comments in read names
    class SeqReader {
    private:
        void inner_read() {
            std::string id, seq;
            std::getline(*stream, id);
            std::getline(*stream, seq);
            if (!id.empty() and !seq.empty())
                next = Contig(Sequence(trim(seq)), trim(id.substr(1, id.size() - 1)));
            else {
                next = Contig();
                return;
            }
            if(fastq) {
                std::getline(*stream, id);
                std::getline(*stream, seq);
            }
        }

        std::istream * stream;
        Contig next;
        bool fastq;
    public:
        explicit SeqReader(const std::string & file_name) {
            if (endsWith(file_name, ".gz")) {
                stream = new gzstream::igzstream(file_name.c_str());
                fastq = endsWith(file_name, "fastq.gz") or endsWith(file_name, "fq.gz");
            } else {
                stream = new std::ifstream(file_name);
                fastq = endsWith(file_name, "fastq") or endsWith(file_name, "fq");
            }
            inner_read();
        }

        Contig read() {
            Contig tmp = std::move(next);
            inner_read();
            return std::move(tmp);
        }

        std::vector<Contig> readAll() {
            std::vector<Contig> res;
            while(!eof()) {
                res.emplace_back(next);
                inner_read();
            }
            return res;
        }

        bool eof() {
            return next.isNull();
        }

        ~SeqReader() {
            free(stream);
        }
    };
}