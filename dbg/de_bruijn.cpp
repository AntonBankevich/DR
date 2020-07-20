//
// Created by anton on 17.07.2020.
//

#pragma once

//
// Created by anton on 08.07.2020.
//
#include <iostream>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>
#include <common/dir_utils.hpp>
#include <queue>
#include "common/cl_parser.hpp"
#include "aligner.hpp"

typedef unsigned __int128 htype;

htype pow(htype base, size_t p) {
    if (p == 0)
        return 1;
    htype tmp = pow(base, p / 2);
    if (p %2 == 1)
        return base * tmp * tmp;
    else
        return tmp * tmp
}

const size_t k = 7000;
const size_t hbase = 239;
const htype kpow = pow(hbase, k - 1);
w = 12000;

class KWH {
public:
    Sequence *seq;
    size_t pos;
    htype hash;
    KWH(Sequence *_seq, size_t _pos): seq(_seq), pos(_pos) {
        hash = 0;
        for(size_t i = _pos; i < _pos + k; i++) {
            hash = hash * hbase + seq->operator[](i);
        }
    }

    KWH(Sequence *_seq, size_t _pos, size_t _hash): seq(_seq), pos(_pos), hash(_hash) {
    }

    bool isMovable() {
        return seq == nullptr;
    }

    KWH next() const {
        htype other_hash = (hash - kpow * dignucl(seq->operator[](pos))) * hbase + seq->operator[](pos + k);
        return KWH(seq, pos + 1, other_hash);
    }

    htype operator<<(char c) const {
        return (hash - kpow * dignucl(seq->operator[](pos))) * hbase + c;
    }
};

class MinQueue {
    std::deque<KWH> q;
public:
    MinQueue() = default;

    void push(const KWH &kwh) {
        while(!q.empty() && q.back().hash >= kwh.hash) {
            q.pop_back();
        }
        q.push_back(kwh);
    }

    void pop(size_t pos) {
        if(!q.empty() && q.front().pos <= pos) {
            q.pop_front();
        }
    }

    bool empty() const {
        return q.empty();
    }

    KWH get() const {
        return q.front();
    }
};

int main(int argc, char **argv) {
    CLParser parser({"reads=", "output-dir=", "threads=8"}, {"o=output-dir", "t=threads"});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    std::cout << "Reading reads" <<std::endl;
    io::SeqReader reader(parser.getValue("reads"));
    std::vector<htype> hashs;
    while(not reader.eof()) {
        Contig read = reader.read();

    }
    std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    cout << std::experimental::filesystem::absolute(dir) << endl;
    std::ofstream os;
    os.open(dir / "log.info");
    std::cout << "Printing results" <<std::endl;
    os.close();
    std::cout << "Finished" <<std::endl;
    return 0;
}