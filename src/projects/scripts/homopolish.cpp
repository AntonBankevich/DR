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

struct DinucleotideInfo {
    size_t len;
//sorted from start to end;
    vector <dinucleotide> dinucleotide_coords;
    map <size_t, size_t> num_dinucleotide;
    vector<vector<size_t>> dinucleotide_read_voting;
    vector<vector<string>> dinucleotide_strings;
    DinucleotideInfo(string & sequence) {

        dinucleotide_coords = GetDinucleotideRepeats(sequence);
        len = sequence.length();
        size_t count = 0;
        for (auto di: dinucleotide_coords ) {
            for (size_t i = 0; i < 2* di.multiplicity; i++)
                num_dinucleotide[di.start + i] = count;
            count++;
            dinucleotide_read_voting.push_back(vector<size_t>());
            dinucleotide_strings.push_back(vector<string>());
        }
    }
    DinucleotideInfo(){
        len = 0;
    }

    string  GenerateDinucleotideConsensus(Logger & logger, string &contig) {
        stringstream res;
        size_t prev_coord = 0;
        size_t total_len = 0;
        for (size_t i = 0; i < dinucleotide_coords.size(); i++) {
            int actual_multiplicity = dinucleotide_coords[i].multiplicity;
            int res_mult = GenerateMajority(i, logger, actual_multiplicity);
            int difference = res_mult - actual_multiplicity;
            total_len += actual_multiplicity * 2;
            res << contig.substr(prev_coord, dinucleotide_coords[i].start - prev_coord);
            prev_coord  = dinucleotide_coords[i].start + 2* actual_multiplicity;
            if (abs(difference) > 0 || actual_multiplicity > 100) {
                logger.info() << "Dimer multiplicity for  pos " << dinucleotide_coords[i].start << endl;
                logger.info() << "Was " << actual_multiplicity << " real " << res_mult << endl;
                stringstream ss;
                for (auto di_mult: dinucleotide_read_voting[i])
                    ss << di_mult << " ";
                logger.info() << ss.str() << endl;
                size_t count = 0;
                for (auto s : dinucleotide_strings[i]) {
                    logger.info() << ">" << count << endl;
                    logger.info() << s << endl;
                    count ++;
                }

            }
            if (abs(difference) > 0 && abs(difference) < dinucleotide::MIN_DINUCLEOTIDE_REPEAT) {
//TODO:: here some shifts;

                actual_multiplicity = res_mult;
            }
            for (size_t j = 0; j < actual_multiplicity; j++) {
                res << dinucleotide_coords[i].seq;
            }
            logger.info() << dinucleotide_coords[i].start << " " << res.str().length() << endl;
        }
        res << contig.substr(prev_coord, len - prev_coord);
        logger.info() << " Totatlly in dimeric regions " << total_len << endl;
        return res.str();
    }
    double GenerateConsensus(size_t i, Logger& logger) {

        size_t total_multiplicity = 0;
        for (size_t mult: dinucleotide_read_voting[i])
            total_multiplicity += mult;
        double avg_mult = (total_multiplicity * 1.0) / dinucleotide_read_voting[i].size();
        return avg_mult;
    }
    size_t GenerateMajority(size_t i, Logger& logger, int actual_multiplicity) {
        sort(dinucleotide_read_voting[i].begin(), dinucleotide_read_voting[i].end());
        size_t max_l = 0;
        size_t res_len = dinucleotide_read_voting[i][0];

        for(auto it = dinucleotide_read_voting[i].begin(); it != dinucleotide_read_voting[i].end();)
        {
            auto r = std::equal_range(it, dinucleotide_read_voting[i].end(), *it);
            int len = std::distance(r.first, r.second);
//Second condition to skip "splitted" dimer alignments that may become majority on laaarge multiplicity dimers
            if ((len >= max_l) && (actual_multiplicity  < dinucleotide::MIN_DINUCLEOTIDE_REPEAT + *it))
            {
                max_l = len;
                res_len = *it;
            }
            it = r.second;
        }
        return res_len;
    }


//alignment for debug only
    bool AddDinucl(size_t pos_center, dinucleotide& di, Logger & logger, AlignmentInfo &aln, string &seq) {
        if (num_dinucleotide.find(pos_center) == num_dinucleotide.end()) {
            logger.info() << "Skipping dinucleotide for contig " << aln.contig_id << " read " << aln.read_id << endl
                          << pos_center << endl;

            return false;
        } else {
            size_t cur_dinucleotide = num_dinucleotide[pos_center];
            logger.info() << " cur dinucleo id" <<  cur_dinucleotide << endl;
            dinucleotide_read_voting[cur_dinucleotide].push_back(di.multiplicity);
            dinucleotide_strings[cur_dinucleotide].push_back(seq);
            return true;
        }
    }

};

struct ContigInfo {
    string sequence;
    string name;
    size_t len;
    vector<uint8_t > quantity;
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
    DinucleotideInfo dinucl;

/*    vector <dinucleotide> dinucleotide_coords;
    vector <int> num_dinucleotide;
    vector<vector<size_t>> dinucleotide_read_voting;
*/
    void FillComplex() {
        size_t current_id = 0;
        size_t total_d = dinucl.dinucleotide_coords.size();
        while (current_id < total_d) {
            size_t cur_start = get_start(dinucl.dinucleotide_coords[current_id]);
            size_t cur_finish = get_finish(dinucl.dinucleotide_coords[current_id]);
            size_t next_id = current_id + 1;
            while ((next_id) < total_d && get_start(dinucl.dinucleotide_coords[next_id]) <= cur_finish) {
                cur_finish = std::max(cur_finish, get_finish(dinucl.dinucleotide_coords[next_id]));
                next_id ++;

            }
            if (cur_finish > len)
                cur_finish = len;
            complex_regions[cur_start] = cur_finish - cur_start;
            complex_strings[cur_start] = vector<string>();
            current_id = next_id;
        }

    }

    ContigInfo(string sequence, string name, string debug_f):sequence(sequence), name(name), debug_f(debug_f), dinucl(sequence) {
        len = sequence.length();
        quantity.resize(len);
        sum.resize(len);
        amounts.resize(len);
        for (size_t i = 0; i < len; i ++){
            sum[i] = 0;
            quantity[i] = 0;
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
                    sort(amounts[i].begin(), amounts[i].end());
                    real_cov = amounts[i][amounts[i].size()/2];
                    cov = int(round(sum[i] * 1.0 / quantity[i]));
                    if (real_cov != cov) {
                        logger.info() << "Median disagree with average at position " << i<<" med/avg: "<< real_cov << "/" << cov << endl;
                        stringstream as;
                        for (auto cc:amounts[i])
                            as << int(cc) << " ";
                        logger.info() << as.str() << endl;
                        cov = real_cov;
                    }
                }

                quantities[quantity[i]]++;
                for (size_t j = 0; j < cov; j++) {
                    ss << sequence[i];
                    debug << total_count << " " << sum[i] << " " << int(quantity[i]) << " " << cov << " " << sequence[i];

                    debug << "    ";
                    for (auto am : amounts[i]){
                        debug << int(am) << " ";
                    }

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
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;
//at least 2



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
            for (auto di: contigs[contig.id].dinucl.dinucleotide_coords){
                logger.info() << di.str() << endl;
            }
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

    void processDinucleotideReadPair (StringContig& read, AlignmentInfo& aln) {
        logger.info() << read.id << endl;
        Sequence uncompressed_read_seq = read.makeSequence();
        read.compress();
        Sequence read_seq = read.makeSequence();
        size_t rlen = read.size();
        if (aln.rc) {
            read_seq = !read_seq;
            uncompressed_read_seq = !uncompressed_read_seq;
            RC(aln);
        }
        logger.info() << aln.alignment_start << " " << aln.alignment_end << endl;
        vector<dinucleotide> dinucleo_read = GetDinucleotideRepeats(read_seq.str());
        size_t used = 0;
        for (auto di: dinucleo_read) {
            size_t pos_start = aln.alignment_start + di.start;
            size_t pos_center = pos_start + di.multiplicity;
            logger.info() << pos_center << " " << di.multiplicity <<  endl;
            logger.info() << di.str()  << endl;
            string uncompressed_dimer = uncompress(di.start, di.start + 2* di.multiplicity, uncompressed_read_seq.str());
            logger.info() << uncompressed_dimer << endl;
            if (contigs[aln.contig_id].dinucl.AddDinucl(pos_center, di, logger, aln, uncompressed_dimer)) {
                used ++;
            }
        }
        if (dinucleo_read.size() > 0) {
            logger.info() << " Used dinucleo " << used << " of " << dinucleo_read.size() << endl;
        }

    }
    bool verifyCigar(vector<cigar_pair> &cigars) {
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
        logger.info() << read.id << endl;
        Sequence read_seq = read.makeSequence();
        size_t rlen = read.size();
        if (aln.rc) {
            read_seq = !read_seq;
            RC(aln);
//            read = read.RC();
        }
//
        string compressed_read = read_seq.str();
        compressed_read.erase(std::unique(compressed_read.begin(), compressed_read.end()), compressed_read.end());
        string contig_seq = contigs[aln.contig_id].sequence.substr(aln.alignment_start , aln.alignment_end - aln.alignment_start);

//        logger.info() << "contig length " << contig_seq.length() << " read length: " << compressed_read.length() << endl;

//        logger.trace() << contig_seq << endl;
//        logger.trace() << compressed_read << endl;


//map-hifi
//asm10
        ;

//strings, match, mismatch, gap_open, gap_extend, width
        auto cigars = align_example( contig_seq.c_str(), compressed_read.c_str(), 1, -5, 5, 2, 50);
        logger.info() << str(cigars) << endl;
        auto str_cigars = str(cigars);
        /*if (!verifyCigar(cigars) || cigars.size() > 100) {
            logger.info() << "STATS: aln length " << aln.length() <<" read length "<< compressed_read.length() << " matched length " << matchedLength(cigars)<<endl;
        } */
        if (matchedLength(cigars) + 300 < compressed_read.length() || !verifyCigar(cigars) || cigars.size() > 300) {
/*            Contig cref(contig_seq, "ref");
            std::vector<Contig > ref = {cref};
            Contig cread(compressed_read, "read");
            RawAligner<Contig> minimapaligner(ref, 1, "ava-hifi");
            auto minimapaln = minimapaligner.align(cread);
            logger.info() <<"minimapped\n";
            logger.info() << minimapaln[0].cigarString() << endl; */
            if (matchedLength(cigars) % 200 == 0) {
                auto full_cigars = align_example( contig_seq.c_str(), compressed_read.c_str(), 1, -5, 5, 2, 500);
                logger.info() << "Check for read " << read.id  << endl;
                logger.info() << "Bandwindth 500\n";
                logger.info() << str(full_cigars) << endl;
                logger.info() << "aln length " << aln.length() <<" read length "<< compressed_read.length() << " matched length " << matchedLength(cigars)<<endl;
            } else {
                logger.info()<< "something wrong but skipping for speedup\n";
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
                for (size_t i = 0; i < (*it).length; i++) {
                    size_t coord = cont_coords + aln.alignment_start + i;
//                    logger.info() << contigs[aln.contig_id].sequence[coord]<< compressed_read[read_coords + i] << endl;
                    if (contigs[aln.contig_id].sequence[coord] == nucl(compressed_read[read_coords + i]) && contigs[aln.contig_id].quantity[coord] != 255 && contigs[aln.contig_id].sum[coord] < 60000) {
                        contigs[aln.contig_id].quantity[coord] ++;
                        contigs[aln.contig_id].sum[coord] += quantities[read_coords + i];
                        if (contigs[aln.contig_id].amounts[coord].size() < 21)
                            contigs[aln.contig_id].amounts[coord].push_back( quantities[read_coords + i]);

                        matches ++;
                    } else {
                        mismatches ++;
                    }
//Only complete traversion of complex regions taken in account;
                    if (contigs[aln.contig_id].complex_regions.find(coord) != contigs[aln.contig_id].complex_regions.end()) {
                        size_t complex_len = contigs[aln.contig_id].complex_regions[coord];
                        if (read_coords + complex_len < matchedLength(cigars)) { //minimapaln[0].seg_from.right) {
                            complex_start = read_coords + i;
                            contig_finish = coord + complex_len;
                            complex_id = coord;
                        }
                    }
                    if (coord == contig_finish) {
                        contigs[aln.contig_id].complex_strings[complex_id].push_back(uncompress(complex_start, read_coords + i, read_seq.str()));
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
//            processDinucleotideReadPair(cur, cur_align);
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


//Awfull copypaste, temporary solution
    void processDinucleo() {

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
//            processReadPair(cur, cur_align);
            processDinucleotideReadPair(cur, cur_align);
//TODO:: appropriate logic for multiple alignment
            do {
                cur_align = readAlignment(compressed_reads);
                aln_count ++;
                if (aln_count % 1000 == 0) {
                    logger.info() << "Processed " << aln_count << " compressed reads " << endl;

                }
                if (cur_compressed == cur_align.read_id) {
//                    processReadPair(cur, cur_align);
                    processDinucleotideReadPair(cur,cur_align);
                }
            } while (cur_compressed == cur_align.read_id);
/*            if (aln_count %2000 == 1999)
                break; */
            cur_compressed = cur_align.read_id;
        }
        logger.info() << "Reads processed, used " <<used_pairs << " filtered by compressed length " << filtered_pairs << endl;
        for (auto& contig: contigs){
            logger.info() << "Generating consensus for contig " << contig.first << endl;
            //contig.second.dinucl.GenerateDinucleotideConsensus(logger);
            // corrected_contigs << ">" << contig.first << '\n' << contig.second.GenerateConsensus(logger) << '\n';
            corrected_contigs << ">" << contig.first << '\n' << contig.second.dinucl.GenerateDinucleotideConsensus(logger, contig.second.sequence) << '\n';
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
