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
#include <spoa/spoa.hpp>
#include "minimap2/ksw2.h"

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


struct dinucleotide {
    size_t start, multiplicity;
    //do we need sequence?
    string seq;
    static const size_t MIN_DINUCLEOTIDE_REPEAT = 10;

    dinucleotide(size_t start, size_t multiplicity, string seq): start(start), multiplicity(multiplicity), seq(seq) {
    }
    string str() {
        stringstream ss;
        ss << start << " " << multiplicity << " " <<seq;
        return ss.str();
    }
};

template <class Contig>
vector<dinucleotide> GetDinucleotideRepeats(const Contig& compressed_contig) {
    size_t len = compressed_contig.size();
    size_t start = 0;
    vector<dinucleotide> res;
    while (start < len) {
        size_t multiplicity = 1;
        while ((start + 2 * multiplicity + 1 < len) &&
               compressed_contig[start] == compressed_contig[start + 2 * multiplicity ] &&
               compressed_contig[start + 1] == compressed_contig[start + 2 * multiplicity + 1]) {
            multiplicity++;
        }

        if (multiplicity >= dinucleotide::MIN_DINUCLEOTIDE_REPEAT)
        {
            stringstream ss;
            ss << compressed_contig[start] << compressed_contig[start+1];
            res.emplace_back(start, multiplicity, ss.str());
//Dirty hack to avoid reporting overlapping dinucleotide repeats
            start ++;
        }
        start = start + 2 * multiplicity - 1;
    }
    return res;
}

struct ContigInfo {
    string sequence;
    string name;
    size_t len;
    vector<uint8_t > quantity;
//array size. possibly switch to []
    static const size_t VOTES_STORED = 21;
    //To avoid multiple resizes real length stored in first element of each vector
    vector<vector<uint8_t>> amounts;
    vector<uint16_t> sum;
    size_t zero_covered = 0;
    string debug_f;
//neighbourhoud for complex regions;
    static const size_t COMPLEX_EPS = 5;
    static const size_t MAX_CONSENSUS_COVERAGE = 25;
    //Complex_regions: map from start to length;
    map<size_t, size_t> complex_regions;
    map<size_t, vector<string>> complex_strings;
//TODO more efficient data structure

/*    vector <dinucleotide> dinucleotide_coords;
    vector <int> num_dinucleotide;
    vector<vector<size_t>> dinucleotide_read_voting;
*/
    void FillComplex() {
        size_t current_id = 0;
        auto dinucleotide_coords = GetDinucleotideRepeats(sequence);
        size_t total_d = dinucleotide_coords.size();
        while (current_id < total_d) {
            size_t cur_start = get_start(dinucleotide_coords[current_id]);
            size_t cur_finish = get_finish(dinucleotide_coords[current_id]);
            size_t next_id = current_id + 1;
            while ((next_id) < total_d && get_start(dinucleotide_coords[next_id]) <= cur_finish) {
                cur_finish = std::max(cur_finish, get_finish(dinucleotide_coords[next_id]));
                next_id ++;

            }
            if (cur_finish > len)
                cur_finish = len;
            complex_regions[cur_start] = cur_finish - cur_start;
            complex_strings[cur_start] = vector<string>();
            current_id = next_id;
        }

    }

    ContigInfo(string sequence, string name, string debug_f):sequence(sequence), name(name), debug_f(debug_f) {
        len = sequence.length();
        quantity.resize(len);
        sum.resize(len);
        amounts.resize(len);
        for (size_t i = 0; i < len; i ++){
            sum[i] = 0;
            quantity[i] = 0;
            amounts[i].resize(VOTES_STORED + 1);
            amounts[i][0] = 0;
        }
        FillComplex();
    }
    ContigInfo() = default;

    size_t get_finish(dinucleotide d) {
        return d.start + COMPLEX_EPS + d.multiplicity * 2;
    }
    size_t get_start(dinucleotide d) {
        if (d.start < COMPLEX_EPS) return 0;
        else return d.start - COMPLEX_EPS;
    }


    void AddRead() {

    }
/*
 * string getSequence(AlignmentInfo al){

    } */

    bool complexStartPosition(size_t ind){
        return (complex_regions.find(ind) != complex_regions.end());
    }

    string MSAConsensus(vector<string> &s) {

//Magic consts from spoa default settings
        auto alignment_engine = spoa::AlignmentEngine::Create(
// -8 in default for third parameter(gap) opening, -6 for forth(gap extension)
                spoa::AlignmentType::kNW, 5, -4, -8, -6, -10, -4);  // linear gaps

        spoa::Graph graph{};
        size_t cov = 0;
        for (const auto& it : s) {
            std::int32_t score = 0;
            auto alignment = alignment_engine->Align(it, graph, &score);
            graph.AddAlignment(alignment, it);
            cov ++;
            if (cov == MAX_CONSENSUS_COVERAGE) {
                break;
            }
        }

        return graph.GenerateConsensus();
    }

    string checkMSAConsensus(string s) {
        s.erase(std::unique(s.begin(), s.end()), s.end());
        if (s.length() <= 2 * COMPLEX_EPS + 2) return "TOO SHORT";
        size_t len = s.length();
        for (size_t i = COMPLEX_EPS + 2; i < len - COMPLEX_EPS; i++) {
            if (s[i] != s[COMPLEX_EPS] && s[i] != s[COMPLEX_EPS + 1])
                return "On position " + to_string(i) + " of " + to_string(len) + " not dimeric nucleo";
        }
        return "";
    }


    string GenerateConsensus(Logger & logger){
        io:stringstream ss;
        vector<int> quantities(256);
        ofstream debug;
        debug.open(debug_f + name);
        debug << name << endl;
        size_t total_count = 0 ;
        for (size_t i = 0; i < len; ) {
            if (complexStartPosition(i))
            {
                auto consensus = MSAConsensus(complex_strings[i]);
                auto check = checkMSAConsensus(consensus);
                logger.info() << "consensus of " << complex_strings[i].size() << ": " << consensus.length()
                << endl << "At position " << ss.str().length() << endl;
                if (check != ""){
                    logger.info() << "SUSPICIOUS CONS " << check << endl;
                    for (auto s: complex_strings[i]) {
                        logger.info() << s << endl;
                    }
                }
                logger.info() << consensus << endl;
                ss << consensus;

                for (size_t j = 0; j < consensus.length(); j++) {
                    debug << total_count << " 0 0 0 " << consensus[j] << endl;
                    total_count ++;
                }

                i += complex_regions[i];
            } else {
                int real_cov = 1;
                int cov = 1;
                if (quantity[i] == 0) {
                    zero_covered++;
                } else {
                    size_t real_len =  amounts[i][0];

/*
                    stringstream uns;
                    uns << " unsorted";
                    for (auto j = 0; j < amounts[i].size(); j++)
                        uns << int(amounts[i][j]) << " ";
                    logger.info() << uns.str() << endl;
*/
                    sort(amounts[i].begin() + 1, amounts[i].begin() + 1 + real_len);
                    real_cov = amounts[i][(1 + real_len)/2];

                    cov = int(round(sum[i] * 1.0 / quantity[i]));
                    if (real_cov != cov) {
                        logger.info() << "Median disagree with average at position " << i<<" med/avg: "<< real_cov << "/" << cov << endl;
                        stringstream as;
                        for (auto j = 0; j < amounts[i][0]; j++)
                            as << int(amounts[i][j + 1]) << " ";
                        logger.info() << as.str() << endl;
                        cov = real_cov;
                    }
                }

                quantities[quantity[i]]++;
                for (size_t j = 0; j < cov; j++) {
                    ss << sequence[i];
                    debug << total_count << " " << sum[i] << " " << int(quantity[i]) << " " << cov << " " << sequence[i];

                    debug << "    ";
                    for (auto jj = 0; jj < amounts[i][0]; jj++)
                        debug << int(amounts[i][jj + 1]) << " ";
                    debug << endl;
                    total_count++;
                }
                i++;
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
    static const size_t SW_BANDWIDTH = 10;
//Used if we see that something wrong is with first
    static const size_t SW_SECOND_BANDWIDTH = 50;
//We do not believe matches on the ends of match region
    static const size_t MATCH_EPS = 1;

    AssemblyInfo (CLParser& parser):parser(parser), logger(), used_pairs(0), filtered_pairs(0), get_uncovered_only(false) {
//TODO paths;
        logging::LoggerStorage ls("/home/lab44/work/homo_compress/", "homopolish");
        logger.addLogFile(ls.newLoggerFile());
        logger.info() <<"reading contigs";
        std::vector<Contig> assembly = io::SeqReader(parser.getValue("contigs")).readAllContigs();
        string debug_f = parser.getValue("debug");
        for (auto contig: assembly) {
//TODO switch to Contig()
            contigs.emplace(contig.id, ContigInfo(contig.seq.str(), contig.id, debug_f));

            logger.info() << contig.id << endl;

        }

    }
    void AddRead(string contig_name){
        contigs[contig_name].AddRead();
    }


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
    size_t getUncompressedPos(const string &s, size_t pos){
        size_t real_pos = 0;
        while (pos > 0) {
            while (s[real_pos] == s[real_pos + 1]) {
                real_pos ++;
            }
            pos --;
            real_pos ++;
        }
        return real_pos;
    }

    string uncompress (size_t from, size_t to, const string &s) {
//        logger.info() << from << " " << to << " " << s.length() << endl;
        size_t real_from = getUncompressedPos(s, from);
//can be optimised
        size_t real_to = getUncompressedPos(s, to);
        if (real_to >= real_from)
            return s.substr(real_from,  real_to -real_from);
        else {
            logger.info() << "BULLSHIT";
            return "";
        }
    }

//Some strange cigars are output
    bool verifyCigar(vector<cigar_pair> &cigars) {
//TODO constant??
        if (cigars.size() > 200) {
            return false;
        }
        for(size_t i = 0; i + 1 < cigars.size(); i++) {
            if (cigars[i].type != 'M' && cigars[i+1].type != 'M') {
                return false;
            }
        }
        return true;
    }

    string str(vector<cigar_pair> &cigars) {
        stringstream ss;
        for(size_t i = 0; i < cigars.size(); i++) {
            ss << cigars[i].length << cigars[i].type;
        }
        return ss.str();
    }

    size_t matchedLength(vector<cigar_pair> &cigars) {
        size_t res = 0;
        for(size_t i = 0; i < cigars.size(); i++) {
            if (cigars[i].type == 'M')
                res += cigars[i].length;
        }
        return res;
    }
    void processReadPair (StringContig& read, AlignmentInfo& aln) {
//        logger.info() << read.id << endl;
        Sequence read_seq = read.makeSequence();
        size_t rlen = read.size();
        if (aln.rc) {
            read_seq = !read_seq;
            RC(aln);
//            read = read.RC();
        }
        ContigInfo& current_contig = contigs[aln.contig_id];
        string compressed_read = read_seq.str();
        compressed_read.erase(std::unique(compressed_read.begin(), compressed_read.end()), compressed_read.end());
        string contig_seq = current_contig.sequence.substr(aln.alignment_start , aln.alignment_end - aln.alignment_start);

//strings, match, mismatch, gap_open, gap_extend, width
        auto cigars = align_example( contig_seq.c_str(), compressed_read.c_str(), 1, -5, 5, 2, SW_BANDWIDTH);

        auto str_cigars = str(cigars);
        /*if (!verifyCigar(cigars) || cigars.size() > 100) {
            logger.info() << "STATS: aln length " << aln.length() <<" read length "<< compressed_read.length() << " matched length " << matchedLength(cigars)<<endl;
        } */
        size_t matched_l = matchedLength(cigars);
        if (matched_l + 300 < aln.length() || !verifyCigar(cigars)) {
/*            Contig cref(contig_seq, "ref");
            std::vector<Contig > ref = {cref};
            Contig cread(compressed_read, "read");
            RawAligner<Contig> minimapaligner(ref, 1, "ava-hifi");
            auto minimapaln = minimapaligner.align(cread);
            logger.info() <<"minimapped\n";
            logger.info() << minimapaln[0].cigarString() << endl; */

//Likely long prefix clipped

            logger.info() << read.id << endl;
            logger.info() << str(cigars) << endl;
            logger.info() << "aln length " << aln.length() <<" read length "<< compressed_read.length() << " matched length " << matched_l<<endl;
            //TODO: we can map to reverse complement, but do not want
            if (matched_l < 100) {
                logger.info() << "Ultrashort alignmnent, doing nothing" << endl;
            } else {
                cigars = align_example( contig_seq.c_str(), compressed_read.c_str(), 1, -5, 5, 2, SW_SECOND_BANDWIDTH);
                logger.info() << "Replaced using bandwindth "<< SW_SECOND_BANDWIDTH << endl;
                logger.info() << str(cigars) << endl;
            }
        }

//        logger.info() << "Aligned " << minimapaln.size()<<"\n";
//        if (minimapaln.size() == 0) {
//            return;
//        }
//        logger.info() << "Shift_to " << minimapaln[0].seg_to.left << " shift from " << minimapaln[0].seg_from.left << " rc" << aln.rc << " " << minimapaln[0].cigarString()<< endl;

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
        size_t cont_coords = 0; //minimapaln[0].seg_to.left;
        size_t read_coords = 0; //minimapaln[0].seg_from.left;
        size_t matches = 0;
        size_t mismatches = 0;
        size_t indels = 0;
        size_t complex_start = -1;
        size_t contig_finish = -1;
        size_t complex_id = -1;
//        logger.info() << "quantities_calc\n";
//        for (auto it = minimapaln[0].begin(); it != minimapaln[0].end(); ++it) {
        for (auto it = cigars.begin(); it != cigars.end(); ++it) {
            if ((*it).type == 'I') {
                read_coords += (*it).length;
                indels += (*it).length;
            } else if ((*it).type == 'D') {
                cont_coords += (*it).length;
//TODO do not forget add complex regions here
                indels += (*it).length;
            } else {
                for (size_t i = MATCH_EPS; i + MATCH_EPS< (*it).length; i++) {
                    size_t coord = cont_coords + aln.alignment_start + i;

//                    logger.info() << current_contig.sequence[coord]<< compressed_read[read_coords + i] << endl;
                    if (current_contig.sequence[coord] == nucl(compressed_read[read_coords + i]) && current_contig.quantity[coord] != 255 && current_contig.sum[coord] < 60000) {
                        current_contig.quantity[coord] ++;
                        current_contig.sum[coord] += quantities[read_coords + i];
                        if (current_contig.amounts[coord][0] < ContigInfo::VOTES_STORED)
                            current_contig.amounts[coord][++current_contig.amounts[coord][0]] = quantities[read_coords + i];

                        matches ++;
                    } else {
                        mismatches ++;
                    }
//Only complete traversion of complex regions taken in account;
                    if (current_contig.complex_regions.find(coord) != current_contig.complex_regions.end()) {
                        size_t complex_len = current_contig.complex_regions[coord];
                        if (read_coords + complex_len < matchedLength(cigars)) {
                            complex_start = read_coords + i;
                            contig_finish = coord + complex_len;
                            complex_id = coord;
                        }
                    }
                    if (coord == contig_finish) {
                        current_contig.complex_strings[complex_id].push_back(uncompress(complex_start, read_coords + i, read_seq.str()));
                        contig_finish = -1;
                    }
                }
                read_coords += (*it).length;
                cont_coords += (*it).length;
            }
        }

//        logger.info() << "Matches "<< matches << " mismatches " << mismatches <<" indels " << indels << endl;
        return;
//ssw align check;
//       aligner.Align(compressed_read.c_str(), contig_seq.c_str(), contig_seq.length(), filter, &alignment, maskLen);
//        logger.info() << "Cigar " << alignment.cigar_string << endl;

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
//                if (reads_count %100 == 99) {
//                    exit(0);
//                }
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
//                    processDinucleotideReadPair(cur,cur_align);
                }
            } while (cur_compressed == cur_align.read_id);
/*            if (aln_count %2000 == 1999)
                break; */
            cur_compressed = cur_align.read_id;
        }
        logger.info() << "Reads processed, used " <<used_pairs << " filtered by compressed length " << filtered_pairs << endl;
        for (auto& contig: contigs){
            logger.info() << "Generating consensus for contig " << contig.first << endl;

            corrected_contigs << ">" << contig.first << '\n' << contig.second.GenerateConsensus(logger) << '\n';
//            corrected_contigs << ">" << contig.first << '\n' << contig.second.dinucl.GenerateDinucleotideConsensus(logger, contig.second.sequence) << '\n';
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

    cout << " parsed";

     AssemblyInfo assemblyInfo(parser);
//    if (parser.getValue("dinucleo") == "true")
//        assemblyInfo.processDinucleo();
//    else
     assemblyInfo.process();

 }
