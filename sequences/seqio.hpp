#pragma once

#include <string>
#include "../common/string_utils.hpp"
#include "stream.hpp"
#include "contigs.hpp"
#include <iterator>

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
        class Iterator : public std::iterator<std::forward_iterator_tag, Contig, size_t, Contig *, Contig&> {
        private:
            SeqReader &reader;
            bool isend;
        public:
            Iterator(SeqReader &_reader, bool _isend) :reader(_reader), isend(_isend) {
                if(reader.eof()) {
                    isend = true;
                }
            }

            void operator++() {
                reader.inner_read();
                if(reader.eof()) {
                    isend = true;
                }
            }
            const Contig &operator*() const {
                return reader.next;
            }

            bool operator==(const Iterator &other) {
                return isend == other.isend;
            }
            bool operator!=(const Iterator &other) {
                return isend != other.isend;
            }
        };

        class SeqIterator : public std::iterator<std::forward_iterator_tag, Sequence, size_t, Sequence *, Sequence&> {
        private:
            SeqReader &reader;
            bool isend;
        public:
            SeqIterator(SeqReader &_reader, bool _isend) :reader(_reader), isend(_isend) {
                if(reader.eof()) {
                    isend = true;
                }
            }

            void operator++() {
                reader.inner_read();
                if(reader.eof()) {
                    isend = true;
                }
            }
            const Sequence &operator*() const {
                return reader.next.seq;
            }
            bool operator==(const SeqIterator &other) {
                return isend == other.isend;
            }
            bool operator!=(const SeqIterator &other) {
                return isend != other.isend;
            }
        };

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

        Iterator begin() {
            return Iterator(*this, false);
        }

        Iterator end() {
            return Iterator(*this, true);
        }

        SeqIterator seqbegin() {
            return SeqIterator(*this, false);
        }

        SeqIterator seqend() {
            return SeqIterator(*this, true);
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