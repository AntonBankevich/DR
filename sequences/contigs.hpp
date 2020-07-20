#pragma once

#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <cstring>
#include <sstream>
#include <unordered_map>
#include "nucl.hpp"
#include "IntrusiveRefCntPtr.h"
#include "verify.hpp"
#include "sequence.hpp"

using std::string;
using std::string;
namespace basic {
    inline string Reverse(const string &s){
        if (s[0] == '-')
            return s.substr(1);
        else
            return "-" + s;
    }
}

template<class T>
class Segment{
public:
    const size_t left;
    const size_t right;
    const T &contig;
    Segment(const T &contig_, size_t left_, size_t right_) : left(left_), right(right_), contig(contig_){
        VERIFY(0 <= left and left <= right and right <= contig.size())
    }

    size_t dist(const Segment<T> &other) const {
        VERIFY(contig == other.contig)
        if (right <= other.left)
            return other.left - right;
        else if (other.right <= left)
            return left - other.right;
        else
            return 0;
    }

    Segment<T> RC() const {
        return Segment(contig.rc(), contig.size() - right, contig.size() - left);
    }

    bool inter(const Segment &other) const {
        return contig == other.contig and not (right <= other.left or left >= other.right);
    }

    int interSize(const Segment &other) const {
        if (not inter(other))
            return -1;
        else
            return int(std::min(right, other.right) - std::max(left, other.left));
    }
};

template <class T>
inline std::ostream& operator<<(std::ostream& os, const Segment<T>& seg)
{
    os << seg.contig.id << "[" << seg.left << ":";
    if (seg.right > seg.contig.size() * 3 / 4)
        os << seg.contig.size() << "-" << (seg.contig.size() - seg.right);
    else
        os << seg.right;
    os << "]";
    return os;
}

template<class T>
class NamedSequence {
public:
    string id;
    Sequence seq;
//protected:
//    T * _rc;
public:
//    NamedSequence(const Sequence &_seq, string _id, T *_rc) : seq(_seq), id(std::move(_id)), _rc(_rc){
//    }

    NamedSequence(const Sequence &_seq, string _id) : seq(_seq), id(std::move(_id)){
//        _rc = new T(!seq, basic::Reverse(id), static_cast<T*>(this));
    }

    Segment<T> asSegment() const {
        return Segment<T>(*this, 0u, size());
    }

    Segment<T> segment(size_t left, size_t right) const {
        return Segment<T>(*this, left, right);
    }

    Segment<T> suffix(size_t pos) const {
        if (pos < 0)
            pos = size() + pos;
        if (pos < 0)
            pos = 0;
        if (pos > size())
            pos = size();
        return Segment<T>(*this, pos, size());
    }

    Segment<T> prefix(size_t len) const {
        len = min(len, size());
        return Segment<T>(*this, 0, len);
    }

    size_t size() const {
        return seq.size();
    }

//    NamedSequence &rc() const {
//        return *_rc;
//    }

    bool operator==(const NamedSequence &other) const {
        return id == other.id;
    }

    char operator[](size_t ind) const {
        return nucl(seq[ind]);
    }

    string str() const {
        return seq.str();
    }

    bool isNull() const {
        return seq.empty();
    }
};

class Contig: public NamedSequence<Contig> {
public:
    Contig(): NamedSequence(Sequence(), ""){
    }

    Contig(const Sequence &_seq, const string &_id): NamedSequence(_seq, _id) {
    }

//    Contig(const Sequence &_seq, const string &_id, Contig *_rc): NamedSequence(_seq, _id, _rc) {
//    }

    Contig(const string &_seq, const string &_id): NamedSequence(Sequence(_seq), _id) {
    }

//    Contig(const string &_seq, const string &_id, Contig *_rc): NamedSequence(Sequence(_seq), _id, _rc) {
//    }
};

template <class T>
class SequenceCollection {
private:
    std::unordered_map<std::string, T *> items;
public:
    SequenceCollection() {
    }

    SequenceCollection(const std::vector<T*> & sequences) {
        for(T *item: sequences){
            items[item->id] = item;
        }
    }
};