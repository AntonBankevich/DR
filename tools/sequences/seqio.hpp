#pragma once

#include "common/string_utils.hpp"
#include "stream.hpp"
#include "contigs.hpp"
#include <iterator>
#include <string>
#include <functional>
#include <utility>

namespace io {

    template<class Reader>
    class ContigIterator : public std::iterator<std::forward_iterator_tag, Contig, size_t, Contig *, Contig&> {
    private:
        Reader &reader;
        bool isend;
    public:
        ContigIterator(Reader &_reader, bool _isend) :reader(_reader), isend(_isend) {
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

        bool operator==(const ContigIterator &other) const {
            return isend == other.isend;
        }
        bool operator!=(const ContigIterator &other) const {
            return isend != other.isend;
        }
    };

    template<class Reader>
    class SeqIterator : public std::iterator<std::forward_iterator_tag, Sequence, size_t, Sequence *, Sequence&> {
    private:
        Reader &reader;
        bool isend;
    public:
        SeqIterator(Reader &_reader, bool _isend) :reader(_reader), isend(_isend) {
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
        bool operator==(const SeqIterator &other) const {
            return isend == other.isend;
        }
        bool operator!=(const SeqIterator &other) const {
            return isend != other.isend;
        }
    };

    //    TODO: Deal with corrupted files, comments in read names
    class SeqReader {
    private:
        void inner_read() {
            std::string id, seq;
            std::getline(*stream, id);
            std::getline(*stream, seq);
            if (!id.empty() and !seq.empty()) {
                seq = trim(seq);
                transform(seq);
                next = Contig(Sequence(seq), trim(id.substr(1, id.size() - 1)));
            }else {
                next = Contig();
                return;
            }
            if(fastq) {
                std::getline(*stream, id);
                std::getline(*stream, seq);
            }
        }

        std::function<void(std::string &)> transform;
        std::istream * stream;
        Contig next;
        bool fastq;
    public:
        friend class ContigIterator<SeqReader>;
        friend class SeqIterator<SeqReader>;

        SeqReader(const std::string & file_name,
                std::function<void(std::string &)> _transform = [](std::string &s) {}) :
                    transform(std::move(_transform)) {
            if (endsWith(file_name, ".gz")) {
                stream = new gzstream::igzstream(file_name.c_str());
                fastq = endsWith(file_name, "fastq.gz") or endsWith(file_name, "fq.gz");
            } else {
                stream = new std::ifstream(file_name);
                fastq = endsWith(file_name, "fastq") or endsWith(file_name, "fq");
            }
            inner_read();
        }

        SeqReader(SeqReader &&other) = default;

        static SeqReader CompressingReader(const std::string &file_name) {
            return {file_name, [](std::string &s) {compress_inplace(s);}};
        }

        ContigIterator<SeqReader> begin() {
            return {*this, false};
        }

        ContigIterator<SeqReader> end() {
            return {*this, true};
        }

        SeqIterator<SeqReader> seqbegin() {
            return {*this, false};
        }

        SeqIterator<SeqReader> seqend() {
            return {*this, true};
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
            return std::move(res);
        }

        bool eof() {
            return next.isNull();
        }

        ~SeqReader() {
            delete stream;
        }
    };

}