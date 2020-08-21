//
// Created by anton on 08.07.2020.
//

#pragma once

#include <vector>
#include <omp.h>
#include "sequences/contigs.hpp"
#include "minimap2/cigar.hpp"
#include "minimap2/minimap_interface.hpp"
namespace cigar_event_type {
    enum CigarEventType : size_t {
        match = 0, ins = 1, del = 2
    };
}

template<class U, class V>
struct CigarAlignment {
    struct CigarEvent {
        size_t from_pos;
        size_t to_pos;
        size_t type;
        size_t length;

        size_t from_end() const {
            if (type == cigar_event_type::del)
                return from_pos;
            else
                return from_pos + length;
        }

        size_t to_end() const {
            if (type == cigar_event_type::ins)
                return to_pos;
            else
                return to_pos + length;
        }
    };

    class Iterator {
    private:
        CigarAlignment<U, V> &al;
        size_t from_pos;
        size_t to_pos;
    public:
        size_t cigar_pos;
        Iterator(CigarAlignment<U, V> &al, size_t from_pos, size_t to_pos, size_t cigar_pos) :
                    al(al), from_pos(from_pos), to_pos(to_pos), cigar_pos(cigar_pos) {
        }

        bool valid() const {
            return from_pos != size_t (-1) && to_pos != size_t(-1);
        }

        CigarEvent operator*() {
            VERIFY(valid());
            size_t block_len = al.cigar_container->cigar[cigar_pos] >> 4u;
            size_t event = al.cigar_container->cigar[cigar_pos] & 15u;
            return CigarEvent(from_pos, to_pos, event, block_len);
        }

        Iterator & operator++() {
            VERIFY(valid());
            size_t block_len = al.cigar_container->cigar[cigar_pos] >> 4u;
            size_t event = al.cigar_container->cigar[cigar_pos] & 15u;
            if (event == cigar_event_type::ins) {
                from_pos += block_len;
            } else if (event == cigar_event_type::del) {
                from_pos += block_len;
            } else {
                from_pos += block_len;
                to_pos += block_len;
            }
            return *this;
        }

        Iterator & operator--() {
            VERIFY(valid());
            cigar_pos -= 1;
            size_t block_len = al.cigar_container->cigar[cigar_pos] >> 4u;
            size_t event = al.cigar_container->cigar[cigar_pos] & 15u;
            if (event == 1) {
                from_pos -= block_len;
            } else if (event == 2) {
                from_pos -= block_len;
            } else {
                from_pos -= block_len;
                to_pos -= block_len;
            }
            return *this;
        }

        bool operator==(const Iterator &other) const {
            return &al == *other.al && cigar_pos == other.cigar_pos;
        }

        bool operator!=(const Iterator &other) const {
            return !this->operator==(other);
        }
    };

    Segment<U> seg_from;
    Segment<V> seg_to;
    mm_extra_t *cigar_container; //For compatibility with minimap

    CigarAlignment(RawAlignment &&other, const U &contig_from, const V &contig_to):
                        seg_from(contig_from, other.seg_from.left, other.seg_from.right),
                        seg_to(contig_to, other.seg_to.left, other.seg_to.right), cigar_container(other.cigar_container) {
//        cout << "CigarAlignment From moved raw: " << cigar_container << endl;
        other.cigar_container = nullptr;
    }

    CigarAlignment(const Segment<U> &segFrom, const Segment<V> &segTo, mm_extra_t *cigarContainer) :
                        seg_from(segFrom), seg_to(segTo), cigar_container(cigarContainer) {}

    CigarAlignment(const CigarAlignment &other) = delete;

    CigarAlignment(CigarAlignment &&other)  noexcept : seg_from(other.seg_from), seg_to(other.seg_to), cigar_container(other.cigar_container) {
        other.cigar_container = nullptr;
    }

    size_t size() const {
        return std::max(seg_from.size(), seg_to.size());
    }

//    Warning. This method can destroy contents alignment. Use it as transformation.
    static std::vector<CigarAlignment<U, V>> split(CigarAlignment<U, V> && al, size_t break_size, size_t min_size) {
        std::vector<CigarAlignment<U, V>> res;
        Iterator start = al.begin();
        Iterator prev = al.begin();
        Iterator finish = al.end();
        for(Iterator cur = al.begin(); cur != finish; ++cur) {
            CigarEvent event = *cur;
            if(event.type == cigar_event_type::match) {
                prev = cur;
                ++prev;
                if(start == finish) {
                    start = cur;
                }
            } else if (prev != finish) {
                size_t len = 0;
                for(size_t i = event.from_pos; i + 1 < event.from_end(); i++) {
                    if(seg_from.contig()[i] != al.seg_from.contig()[i + 1])
                        len += 1;
                }
                for(size_t i = event.to_pos; i + 1 < event.to_end(); i++) {
                    if(seg_to.contig()[i] != al.seg_to.contig()[i + 1])
                        len += 1;
                }
                if(len > break_size) {
                    CigarAlignment<U, V> next = al.subAlignment(start, prev);
                    if(next.size() >= min_size)
                        res.push_back(al.subAlignment(start, prev));
                    start = finish;
                    prev = finish;
                }
            }
        }
        if(start != finish) {
            if(start == al.begin()) {
                if(al.size() > min_size)
                    res.emplace_back(al);
            } else {
                CigarAlignment<U, V> next = al.subAlignment(start, prev);
                if(next.size() >= min_size)
                    res.push_back(al.subAlignment(start, prev));
            }
        }
        return res;
    }

//    Subalignment is always initially bound at M blocks in cigar.
//    We shrink subalignment further to satisfy the condition that the first and last nucleotides are matching
    CigarAlignment<U, V> subAlignment(Iterator it_left, Iterator it_right) const {
        size_t cut_left = 0;
        size_t cut_right = 0;
        Segment<U> new_from(seg_from.contig(), it_left.from_pos, it_right.from_pos);
        Segment<V> new_to(seg_to.contig(), it_left.to_pos, it_right.to_pos);
        {
            bool found = false;
            while (!found && it_left != it_right) {
                CigarEvent event = *it_left;
                if (event.type == cigar_event_type::match) {
                    for (size_t i = 0; i < event.length; i++) {
                        if (seg_from.contig()[event.from_pos + i] == seg_to.contig()[event.to_pos + i]) {
                            cut_left = i;
                            new_from.left = event.from_pos + i;
                            new_to.left = event.to_pos + i;
                            found = true;
                            break;
                        }
                    }
                }
                ++it_left;
            }
            VERIFY(found);
        }
        {
            bool found = false;
            while(it_left != it_right) {
                --it_right;
                CigarEvent event = *it_right;
                if (event.type == cigar_event_type::match) {
                    for (size_t i = 0; i < event.length; i++) {
                        if (seg_from.contig()[event.from_pos + event.length - 1 - i] ==
                                seg_to.contig()[event.to_pos + event.length - 1 - i]) {
                            new_from.right = event.from_pos + event.length - i;
                            new_to.right = event.to_pos + event.length - i;
                            cut_right = i;
                            found = true;
                            break;
                        }
                    }
                }
            }
            VERIFY(found);
        }
        uint32_t capacity = it_right.cigar_pos - it_left.cigar_pos + 1 + sizeof(mm_extra_t)/4;
        auto *subcigar = (mm_extra_t*)calloc(capacity, 4);
        subcigar->capacity = capacity;
        subcigar->n_cigar = it_right.cigar_pos - it_left.cigar_pos + 1;
        for(size_t i = it_left.cigar_pos; i <= it_right.cigar_pos; i++) {
            subcigar->cigar[i - it_left.cigar_pos] = cigar_container->cigar[i];
        }
        subcigar->cigar[0] -= cut_left >> 4u;
        subcigar->cigar[subcigar->n_cigar - 1] -= cut_right >> 4u;
        return {new_from, new_to, subcigar};
    }

    Iterator begin() const {
        return Iterator(*this, seg_from.left, seg_to.left, 0);
    }

    Iterator end() const {
        return Iterator(*this, seg_from.right, seg_to.right, cigar_container->n_cigar);
    }

    std::string cigarString() const {
        std::stringstream ss;
        for(Iterator it = begin(); it != end(); ++it) {
            CigarEvent event = *it;
            ss << event.length << "MID"[event.type];
        }
        return ss.str();
    }

    ~CigarAlignment() {
        if (cigar_container != nullptr) {
//            cout << "CigarAlignment destructor" << cigar_container << endl;
            free(cigar_container);
//            cout << "CigarAlignment destructor done" << cigar_container << endl;
        }
    }
};

template<class U, class V, class A>
class AlignmentTranslator {
public:
    AlignmentTranslator() = default;

    virtual A translateOne(CigarAlignment<U, V> &&al) const = 0;

    std::vector<A> translate(std::vector<CigarAlignment<U, V>> &als) {
        std::vector<A> res;
        for(CigarAlignment<U, V> &al: als) {
            res.emplace_back(std::move(translateOne(std::move(al))));
        }
        return std::move(res);
    }

    std::vector<std::vector<A>> translate(std::vector<std::vector<CigarAlignment<U, V>>> &als, size_t threads_num = 1) {
        std::vector<std::vector<A>> res;
        res.resize(als.size());
        omp_set_dynamic(0);
        omp_set_num_threads(1);
#pragma omp parallel for default(none) shared(als, res)
        for(size_t i = 0; i < als.size(); i++) {
            res[i] = translate(als[i]);
        }
        return std::move(res);
    }

    virtual ~AlignmentTranslator() = default;
};