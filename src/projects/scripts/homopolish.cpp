#include <alignment/aligner.hpp>
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <unordered_set>
#include <vector>
#include <iostream>

using namespace std;

struct AlignmentInfo {
    string read_id;
    string contig_id;
    size_t alignment_start;
    size_t alignment_end;
    bool rc;
    int length() {
        return alignment_start - alignment_end;
    }
};


struct ContigInfo {
    string sequence;
    string name;
    size_t len;
    vector<int8_t > quantity;
    vector<int16_t> sum;

    ContigInfo(string sequence, string name):sequence(sequence), name(name) {
        len = sequence.length();
        quantity.resize(len);
        sum.resize(len);
        for (size_t i = 0; i < len; i ++){
            sum[i] = 0;
            quantity[i] = 0;
        }
    }
    ContigInfo() = default;

    void AddRead(){

    }
/*
 * string getSequence(AlignmentInfo al){

    } */
    string GenerateConsensus(){
        io:stringstream ss;
        for (size_t i = 0; i < len; i++) {
            int cov = 1;
            if (quantity[i] == 0) {
                cout << "EMPTY COVERAGE" ;
//ss <<"EMPTY"
//              exit(1);
            }

            cov = trunc(sum[i] / quantity[i]);
            for (size_t j = 0; j < cov; j++)
                ss << sequence[i];
        }
        return ss.str();
    }

};
struct AssemblyInfo {
    map<string, ContigInfo> contigs;
    CLParser parser;
    AssemblyInfo (CLParser parser):parser(parser) {
        std::vector<Contig> assembly = io::SeqReader(parser.getValue("contigs")).readAllContigs();
        for (auto contig: assembly) {
//TODO switch to Contig()
            contigs[contig.id] = ContigInfo(contig.id, contig.seq.str());
        }

    }
    void AddRead(string contig_name){
        contigs[contig_name].AddRead();
    }

    string getSequence(AlignmentInfo al) {
        if (!al.rc) {
            return contigs[al.contig_id].sequence.substr(al.alignment_start, al.alignment_end);
        } else {
//TOFILL
            return "";
        }
    }

    AlignmentInfo readAssembly(std::ifstream& ss){
        AlignmentInfo res;
        if (ss.eof())
            return res;
        ss >> res.read_id >> res.contig_id >> res.alignment_start>>res.alignment_end;
        if (res.contig_id[0] == '-') {
    //1:
            res.contig_id = res.contig_id.substr(1);
            res.rc = true;
        } else {
            res.rc = false;
        }

        return res;
    }

    void processReadPair (StringContig& read, AlignmentInfo aln) {
        string seq = read.seq;
        seq.erase(std::unique(seq.begin(), seq.end()), seq.end());
        if (seq.length() == aln.length()) {
//TODO appropriate condition?
            int cur_ind = 0;
            int shift = 0;
            size_t rlen = read.size();
            int next_ind = cur_ind + 1;
            while (cur_ind < rlen) {
                while (next_ind < rlen && read.seq[cur_ind] == read.seq[next_ind]) {
                    next_ind ++;
                }
            }
            int vote = next_ind - cur_ind;
//TODO: overfill
            contigs[aln.contig_id].quantity[shift + aln.alignment_start] ++;
            contigs[aln.contig_id].sum[shift + aln.alignment_start] += vote;
            cur_ind = next_ind;
            next_ind = cur_ind + 1;
            shift++;
        }
    }

    void process() {
        ifstream reads;
        ifstream compressed_reads;
        io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
//        io::Library lib = oneline::initialize<std::experimental::filesystem::path>("reads.fasta");
        io::SeqReader reader(lib);
        AlignmentInfo cur_align = readAssembly(compressed_reads);
        string cur_compressed = cur_align.read_id;
        string cur_read = "";
        string cur_seq = "";
        StringContig cur;
        while (!compressed_reads.eof()) {
            while (cur.id != cur_compressed) {
                cur = reader.read();
            }

            processReadPair(cur, cur_align);
//TODO:: appropriate logic for multiple alignment
            do {
                cur_align = readAssembly(compressed_reads);
            } while (cur_compressed == cur_align.contig_id);
            cur_compressed = cur_align.contig_id;
        }
        for (auto& contig: contigs){

            cout << contig.first << "/n" << contig.second.GenerateConsensus();
        }
    }
};



/*
 * я выдам записи вида read_id contig_id alignment_start alignment_end
При этом прикладывания на другой стренд будут закодированы в contig_id. К нему в начале будет добавлен "-"
и позиции будут как в реверс комплиментарной строке, а не как в исходной
 *
 */
int main(int argc, char **argv) {
    CLParser parser({"aligned=", "contigs="}, {"reads"},
                    {});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }
    AssemblyInfo assemblyInfo(parser);
    assemblyInfo.process();
/*
    std::unordered_map<std::string, Sequence> ref;
    std::ifstream is;
    is.open(parser.getValue("paf"));
    const std::experimental::filesystem::path out(parser.getValue("output"));
    std::ofstream os;
    std::string line;
    std::unordered_set<std::string> read_set;
    bool check_size = parser.getCheck("full");
    while (std::getline(is, line))
    {
        std::vector<std::string> ss = split(line);
        if(ss[5] != parser.getValue("chr"))
            continue;
        if (check_size) {
            size_t from = std::stoull(ss[7]);
            size_t to = std::stoull(ss[8]);
            size_t sz = std::stoull(ss[1]);
            if (to - from > std::max(std::min(sz * 9 / 10, sz - 1000), sz * 2 / 3)) {
                read_set.insert(ss[0]);
            }
        }
    }
    is.close();
    os.open(out);
    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::SeqReader reader(lib);
    for (StringContig read : reader) {
        if (read_set.find(read.id) != read_set.end()) {
            os << ">" << read.id << "\n" << read.seq << "\n";
        }
    }
    os.close();
    return 0;*/
}