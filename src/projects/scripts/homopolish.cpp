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
using logging::Logger;


struct AlignmentInfo {
    string read_id;
    string contig_id;
    size_t alignment_start;
    size_t alignment_end;
    bool rc;
    int length() {
        return alignment_end - alignment_start;
    }
    string str() {
        stringstream ss;
        ss << read_id << " " << contig_id << " " <<alignment_start << " " << alignment_end << " "<< rc;
        return ss.str();
    }
};


struct ContigInfo {
    string sequence;
    string name;
    size_t len;
    vector<uint8_t > quantity;
    vector<uint16_t> sum;

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
    string GenerateConsensus(Logger& logger){
        io:stringstream ss;
        vector<int> quantities(256);
        for (size_t i = 0; i < len; i++) {
            int cov = 1;
            if (quantity[i] == 0) {

            }
            else {
                cov = trunc(sum[i] / quantity[i]);
            }
            quantities[quantity[i]] ++;
            for (size_t j = 0; j < cov; j++)
                ss << sequence[i];
        }
        logger.info() <<"Contig " << name << endl;
//Constants;
        for (size_t i = 0; i < 20; i++) {
            logger.info() << i << " " << quantities[i] << endl;
        }
        return ss.str();
    }

};
struct AssemblyInfo {
    map<string, ContigInfo> contigs;
    CLParser parser;
    Logger logger;
    size_t used_pairs;
    size_t filtered_pairs;
    AssemblyInfo (CLParser& parser):parser(parser), logger(), used_pairs(0), filtered_pairs(0) {
//TODO paths;
        logging::LoggerStorage ls("/home/lab44/work/homo_compress/", "homopolish");
        logger.addLogFile(ls.newLoggerFile());
        std::vector<Contig> assembly = io::SeqReader(parser.getValue("contigs")).readAllContigs();
        for (auto contig: assembly) {
//TODO switch to Contig()
            contigs[contig.id] = ContigInfo(contig.seq.str(), contig.id);
        }

    }
    void AddRead(string contig_name){
        contigs[contig_name].AddRead();
    }

    /*string getSequence(AlignmentInfo al) {
        if (!al.rc) {
            return contigs[al.contig_id].sequence.substr(al.alignment_start, al.alignment_end - al.alignment_end);
        } else {
//TOFILL
            return "";
        }
    } */

    AlignmentInfo readAlignment(std::ifstream &ss){
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

    void processReadPair (StringContig& read, AlignmentInfo& aln) {
        Sequence read_seq = read.makeSequence();
        size_t rlen = read.size();
        if (aln.rc) {
            read_seq = !read_seq;
            size_t contig_len = contigs[aln.contig_id].len;
            size_t tmp = contig_len - aln.alignment_end;
            aln.alignment_end = contig_len - aln.alignment_start;
            aln.alignment_start = tmp;
//
//            read = read.RC();
        }

        string seq = read_seq.str();
        seq.erase(std::unique(seq.begin(), seq.end()), seq.end());
//TODO: rc reads;
        if (seq.length() == aln.length()) {
//TODO appropriate condition?
            string contig_seq = contigs[aln.contig_id].sequence.substr(aln.alignment_start , aln.alignment_end - aln.alignment_start);
            size_t bad_pos = 0;

            for (size_t i = 0; i < contig_seq.length(); i++)
                if (seq[i]!=contig_seq[i])
                    bad_pos ++;
            if (rlen > 100 && bad_pos * 4 > rlen) {
                logger.info()<< "Lots of bad positions in read " << float(bad_pos)/float(rlen) << endl;
                logger.info ()<< aln.str() << endl;
                filtered_pairs ++;
                return;
            }
/*
 *          logger.info() << seq << endl;
            logger.info()<< aln.read_id << " " << aln.contig_id << endl;
            logger.info() << aln.alignment_start << " " << aln.alignment_end << endl;
            logger.info() << contigs[aln.contig_id].sequence.substr(aln.alignment_start , aln.alignment_end - aln.alignment_start) << endl;
/*            exit(0); */
            int cur_ind = 0;
            int shift = 0;
            int next_ind = cur_ind + 1;
            while (cur_ind < rlen) {
                while (next_ind < rlen && read_seq[cur_ind] == read_seq[next_ind]) {
                    next_ind ++;
                }

                int vote = next_ind - cur_ind;
        //TODO: consts, overfill
                size_t coord = shift + aln.alignment_start;
/*                if (contigs[aln.contig_id].sequence[coord] != nucl(read_seq[cur_ind]))
                    bad_pos ++;
                logger.info() << contigs[aln.contig_id].sequence[coord] << nucl(read_seq[cur_ind]) << endl;
*/
                if (contigs[aln.contig_id].quantity[coord] != 255 && contigs[aln.contig_id].sum[coord] < 60000) {
                    contigs[aln.contig_id].quantity[coord] ++;
                    contigs[aln.contig_id].sum[coord] += vote;
                } else {
                    logger.info() << " overfilled " << contigs[aln.contig_id].quantity[coord] << contigs[aln.contig_id].sum[coord];
                }
                cur_ind = next_ind;
                next_ind = cur_ind + 1;
                shift++;
            }
            used_pairs ++;
            if (rlen > 100 && bad_pos * 4 > rlen) {
                logger.info()<< "Lots of bad positions in read " << float(bad_pos)/float(rlen) << endl;
                logger.info ()<< aln.str() << endl;
            }
        } else {
            filtered_pairs ++;
       //     logger.trace() << seq.length() << " " << aln.length () << endl;
        }
    }

    void process() {

        ifstream compressed_reads;
        ofstream corrected_contigs;
        corrected_contigs.open(parser.getValue("output"));
        compressed_reads.open(parser.getValue("aligned"));
        io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
//        io::Library lib = oneline::initialize<std::experimental::filesystem::path>("reads.fasta");
        io::SeqReader reader(lib);
        logger.info() << "Initialized\n";

        AlignmentInfo cur_align = readAlignment(compressed_reads);
        string cur_compressed = cur_align.read_id;
        string cur_read = "";
        string cur_seq = "";
        StringContig cur;
        logger.info() << "Assembly read\n";
        size_t reads_count = 0;
        size_t aln_count = 1;
        while (!compressed_reads.eof()) {
            bool reads_over = false;
            while (cur.id != cur_compressed) {

                if (reader.eof()) {
                    logger.info() << "Reads are over\n";
                    reads_over = true;
                    break;
                }
                cur = reader.read();
                reads_count ++;
                if (reads_count % 1000 == 0) {
                    logger.info() << "Processed " << reads_count << " original reads " << endl;

                }
            }
            if (reads_over) {
                break;
            }
            processReadPair(cur, cur_align);
//TODO:: appropriate logic for multiple alignment
            do {
                cur_align = readAlignment(compressed_reads);
                aln_count ++;
                if (aln_count % 1000 == 0) {
                    logger.info() << "Processed " << aln_count << " compressed reads " << endl;
                }
            } while (cur_compressed == cur_align.read_id);
            cur_compressed = cur_align.read_id;

        }
        logger.info() << "Reads processed, used " <<used_pairs << " filtered by compressed length " << filtered_pairs << endl;
        for (auto& contig: contigs){

            corrected_contigs << ">" << contig.first << '\n' << contig.second.GenerateConsensus(logger) << '\n';
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
    CLParser parser({"aligned=", "contigs=", "output="}, {"reads"},
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