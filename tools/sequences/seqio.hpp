#pragma once

#include "common/string_utils.hpp"
#include "stream.hpp"
#include "contigs.hpp"
#include <experimental/filesystem>
#include <iterator>
#include <string>
#include <utility>
#include <vector>
#include <functional>
#include <utility>

namespace io {


    typedef std::vector<std::experimental::filesystem::path> Library;

    template<class Reader>
    class ContigIterator {
    private:
        Reader &reader;
        bool isend;
    public:
        typedef StringContig value_type;

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
        StringContig operator*() {
            return std::move(reader.next);
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
        const Sequence &operator*() {
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
    public:
        typedef ContigIterator<SeqReader> Iterator;
    private:
        void inner_read() {
            while (stream != nullptr){
                std::string id, seq;
                std::getline(*stream, id);
                std::getline(*stream, seq);
                trim(seq);
                if (!id.empty() and !seq.empty()) {
                    std::stringstream ss;
                    ss << seq;
                    size_t cnt = 0;
                    while(stream->peek() != EOF && stream->peek() != '>' && stream->peek() != '+') {
                        std::getline(*stream, seq);
                        trim(seq);
                        if (seq.empty())
                            break;
                        ss << seq;
                        cnt += 1;
                    }
                    next = {ss.str(), std::move(trim(id.substr(1, id.size() - 1)))};
                    if(fastq) {
                        std::getline(*stream, id);
                        size_t qlen= 0;
                        while(!stream->eof() && qlen < next.size()) {
                            std::getline(*stream, seq);
                            trim(seq);
                            qlen += seq.size();
                            if (seq.empty())
                                break;
                        }
                    }
                    return;
                }
                nextFile();
            }
            next = StringContig();
        }

        void nextFile() {
            delete stream;
            if (file_it == lib.end()) {
                stream = nullptr;
            } else {
                std::experimental::filesystem::path file_name = *file_it;
                if (endsWith(file_name, ".gz")) {
                    stream = new gzstream::igzstream(file_name.c_str());
                    fastq = endsWith(file_name, "fastq.gz") or endsWith(file_name, "fq.gz");
                } else {
                    stream = new std::ifstream(file_name);
                    fastq = endsWith(file_name, "fastq") or endsWith(file_name, "fq");
                }
                ++file_it;
            }
        }

        std::function<void(std::string &)> transform;
        const Library lib;
        Library::const_iterator file_it;
        std::istream * stream{};
        StringContig next{};
        bool fastq{};
    public:
        friend class ContigIterator<SeqReader>;
        friend class SeqIterator<SeqReader>;

        explicit SeqReader(Library _lib, std::function<void(std::string &)> _transform = [](std::string &s) {}) :
                    transform(std::move(_transform)), lib(std::move(_lib)), file_it(lib.begin()) {
            nextFile();
            inner_read();
        }

        void reset() {
            file_it = lib.begin();
            nextFile();
            inner_read();
        }

        explicit SeqReader(const std::experimental::filesystem::path & file_name,
                  std::function<void(std::string &)> _transform = [](std::string &s) {}) :
                  SeqReader(Library({file_name}), std::move(_transform)) {
        }

        SeqReader(SeqReader &&other) = default;

//        static SeqReader CompressingReader(const std::string &file_name) {
//            return SeqReader(file_name, [](std::string &s) {compress_inplace(s);});
//        }
//
//        static SeqReader CompressingReader(const Library &lib_) {
//            return SeqReader(lib_, [](std::string &s) {compress_inplace(s);});
//        }

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

        StringContig read() {
            StringContig tmp = std::move(next);
            inner_read();
            return std::move(tmp);
        }

        std::vector<StringContig> readAll() {
            std::vector<StringContig> res;
            while(!eof()) {
                res.emplace_back(std::move(next));
                inner_read();
            }
            return std::move(res);
        }

        std::vector<Contig> readAllContigs() {
            std::vector<Contig> res;
            while(!eof()) {
                res.emplace_back(std::move(next.makeContig()));
                inner_read();
            }
            return std::move(res);
        }

        std::vector<Contig> readAllCompressedContigs() {
            std::vector<Contig> res;
            while(!eof()) {
                res.emplace_back(std::move(next.makeCompressedContig()));
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