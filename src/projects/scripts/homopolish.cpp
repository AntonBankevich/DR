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
#include <ssw/ssw_cpp.h>

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
    size_t zero_covered = 0;
    string debug_f;

    ContigInfo(string sequence, string name, string debug_f):sequence(sequence), name(name), debug_f(debug_f) {
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
        ofstream debug;
        debug.open(debug_f + name);
        debug << name << endl;
        size_t total_count = 0 ;
        for (size_t i = 0; i < len; i++) {
            int cov = 1;
            if (quantity[i] == 0) {
                zero_covered ++;
            }

            else {
                cov = int (round(sum[i] * 1.0 / quantity[i]));
            }

            quantities[quantity[i]] ++;
            for (size_t j = 0; j < cov; j++) {
                ss << sequence[i];
                debug << total_count << " " << sum[i] << " " << int(quantity[i]) << " " <<cov << " " << sequence[i] << endl;
                total_count ++;
            }
        }
        logger.info() <<" Contig " << name << " " << sequence.length() << endl;
        logger.info() << "Zero covered (after filtering) " << zero_covered << endl;
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
    bool get_uncovered_only;
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;
//at least 2
    const size_t MIN_DINUCLEOTIDE_REPEAT = 10;
    struct dinucleotide {
        size_t start, multiplicity;
        //do we need sequence?
        dinucleotide(size_t start, size_t multiplicity): start(start), multiplicity(multiplicity){
        }
    };


    AssemblyInfo (CLParser& parser):parser(parser), logger(), used_pairs(0), filtered_pairs(0), get_uncovered_only(false) {
//TODO paths;
        logging::LoggerStorage ls("/home/lab44/work/homo_compress/", "homopolish");
        logger.addLogFile(ls.newLoggerFile());
        std::vector<Contig> assembly = io::SeqReader(parser.getValue("contigs")).readAllContigs();
        string debug_f = parser.getValue("debug");
        for (auto contig: assembly) {
//TODO switch to Contig()
            contigs[contig.id] = ContigInfo(contig.seq.str(), contig.id, debug_f);
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

    void RC(AlignmentInfo & aln){
        size_t contig_len = contigs[aln.contig_id].len;
        size_t tmp = contig_len - aln.alignment_end;
        aln.alignment_end = contig_len - aln.alignment_start;
        aln.alignment_start = tmp;
    }

    vector<dinucleotide> GetDinucleotideRepeats(Contig& compressed_contig){
        size_t len = compressed_contig.size();
        size_t start = 0;
        vector<dinucleotide> res;
        while (start < len) {
            size_t multiplicity = 1;
            while ((start + 2 * multiplicity + 3 < len) &&
                    compressed_contig[start] == compressed_contig[start + 2 * multiplicity + 2] &&
                    compressed_contig[start + 1] == compressed_contig[start + 2 * multiplicity + 3]) {
                multiplicity++;
            }
            if (multiplicity >= MIN_DINUCLEOTIDE_REPEAT) {
                res.emplace_back(start, multiplicity);
            }
            start = start + 2 * multiplicity - 1;
        }
        return res;
    }

    void processReadPair (StringContig& read, AlignmentInfo& aln) {
        logger.info() << read.id << endl;
        Sequence read_seq = read.makeSequence();
        size_t rlen = read.size();
        if (aln.rc) {
            read_seq = !read_seq;
            RC(aln);
//            read = read.RC();
        }
        logger.info() << aln.alignment_start << " " << aln.alignment_end << endl;
        string compressed_read = read_seq.str();
        compressed_read.erase(std::unique(compressed_read.begin(), compressed_read.end()), compressed_read.end());
        string contig_seq = contigs[aln.contig_id].sequence.substr(aln.alignment_start , aln.alignment_end - aln.alignment_start);
        size_t maskLen = 15;

        Contig cref(contig_seq, "ref");
        std::vector<Contig > ref = {cref};
//        logger.info() << "contig length " << contig_seq.length() << " read length: " << compressed_read.length() << endl;

//        logger.trace() << contig_seq << endl;
//        logger.trace() << compressed_read << endl;
        Contig cread(compressed_read, "read");
        std::vector<Contig > reads = {cread};
//map-hifi
//asm10
        RawAligner<Contig> minimapaligner(ref, 1, "ava-hifi");
        auto minimapaln = minimapaligner.align(cread);
//        logger.info() << "Aligned " << minimapaln.size()<<"\n";
        if (minimapaln.size() == 0) {
            return;
        }
        logger.info() << "Shift_to " << minimapaln[0].seg_to.left << " shift from " << minimapaln[0].seg_from.left << " " << minimapaln[0].cigarString()<< endl;

        int cur_ind = 0;
        std:vector<size_t> quantities;
        quantities.resize(compressed_read.length());

        for (size_t i = 0; i < compressed_read.size(); i++){
            size_t count = 0;
            while (cur_ind < read_seq.size() && nucl(read_seq[cur_ind]) == compressed_read[i]) {
                count ++;
                cur_ind ++;
            }
            quantities[i] = count;
        }
        size_t cont_coords = minimapaln[0].seg_to.left;
        size_t read_coords = minimapaln[0].seg_from.left;
        size_t matches = 0;
        size_t mismatches = 0;
        size_t indels = 0;
//        logger.info() << "quantities_calc\n";
        for (auto it = minimapaln[0].begin(); it != minimapaln[0].end(); ++it) {
            if ((*it).type == 1) {
                read_coords += (*it).length;
                indels += (*it).length;
            } else if ((*it).type == 2) {
                cont_coords += (*it).length;
                indels += (*it).length;
            } else {
                for (size_t i = 0; i < (*it).length; i++) {
                    size_t coord = cont_coords + aln.alignment_start + i;
//                    logger.info() << contigs[aln.contig_id].sequence[coord]<< compressed_read[read_coords + i] << endl;
                    if (contigs[aln.contig_id].sequence[coord] == nucl(compressed_read[read_coords + i]) && contigs[aln.contig_id].quantity[coord] != 255 && contigs[aln.contig_id].sum[coord] < 60000) {
                        contigs[aln.contig_id].quantity[coord] ++;
                        contigs[aln.contig_id].sum[coord] += quantities[read_coords + i];
                        matches ++;
                    } else {
                        mismatches ++;
                    }
                }
                read_coords += (*it).length;
                cont_coords += (*it).length;
            }
        }

//        logger.info() << "Matches "<< matches << " mismatches " << mismatches <<" indels " << indels << endl;
        return;
//ssw align check;
        aligner.Align(compressed_read.c_str(), contig_seq.c_str(), contig_seq.length(), filter, &alignment, maskLen);
        logger.info() << "Cigar " << alignment.cigar_string << endl;
//old no_align;
        if (compressed_read.length() == aln.length()) {
            size_t bad_pos = 0;

            for (size_t i = 0; i < contig_seq.length(); i++)
                if (compressed_read[i]!=contig_seq[i])
                    bad_pos ++;
            if (rlen > 100 && bad_pos * 4 > rlen) {
//                logger.info()<< "Lots of bad positions in read " << float(bad_pos)/float(rlen) << endl;
//                logger.info ()<< aln.str() << endl;
                filtered_pairs ++;
                return;
            }
/*
 *          logger.info() << compressed_read << endl;
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

                if (contigs[aln.contig_id].sequence[coord] == nucl(read_seq[cur_ind]) && contigs[aln.contig_id].quantity[coord] != 255 && contigs[aln.contig_id].sum[coord] < 60000) {
                    contigs[aln.contig_id].quantity[coord] ++;
                    contigs[aln.contig_id].sum[coord] += vote;
                } else {
//                    logger.info() << " overfilled " << contigs[aln.contig_id].quantity[coord] << contigs[aln.contig_id].sum[coord];
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
       //     logger.trace() << compressed_read.length() << " " << aln.length () << endl;
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
                if (cur_compressed == cur_align.read_id) {
                    processReadPair(cur, cur_align);
                }
            } while (cur_compressed == cur_align.read_id);
/*            if (aln_count %2000 == 1999)
                break; */
            cur_compressed = cur_align.read_id;

        }
        logger.info() << "Reads processed, used " <<used_pairs << " filtered by compressed length " << filtered_pairs << endl;
        for (auto& contig: contigs){

            corrected_contigs << ">" << contig.first << '\n' << contig.second.GenerateConsensus(logger) << '\n';
        }
        size_t total_zero_covered = 0;
        for (auto & contig: contigs) {
            total_zero_covered += contig.second.zero_covered;
        }
        logger.info() << "Total zero covered "  << total_zero_covered << endl;
    }
};



/*
 * я выдам записи вида read_id contig_id alignment_start alignment_end
При этом прикладывания на другой стренд будут закодированы в contig_id. К нему в начале будет добавлен "-"
и позиции будут как в реверс комплиментарной строке, а не как в исходной
 *
 */
int main(int argc, char **argv) {
    CLParser parser({"aligned=", "contigs=", "output=", "debug="}, {"reads"},
                    {});
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }

    AssemblyInfo assemblyInfo(parser);
    assemblyInfo.process();
}