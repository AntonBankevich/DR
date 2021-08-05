#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <unordered_set>
#include <vector>
#include <iostream>

#include <spoa/spoa.hpp>
#include <ksw2/ksw_wrapper.hpp>


using namespace std;

using logging::Logger;

static const size_t COMPLEX_EPS = 5;
static const size_t MIN_DISTANCE_COMPLEX_REGIONS = 10;
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

//complex region positions sorted
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

size_t get_finish(dinucleotide d) {
    return d.start + COMPLEX_EPS + d.multiplicity * 2;
}

size_t get_start(dinucleotide d) {
    if (d.start < COMPLEX_EPS) return 0;
    else return d.start - COMPLEX_EPS;
}


//extend a bit and merge intersecting regions
template <class Contig>
vector<pair<size_t, size_t>> GetComplexRegions(const Contig& compressed_contig) {
    size_t current_id = 0;
    size_t len = compressed_contig.length();
    vector<pair<size_t, size_t>> complex_regions;
    auto dinucleotide_coords = GetDinucleotideRepeats(compressed_contig);
    size_t total_d = dinucleotide_coords.size();
    while (current_id < total_d) {
        size_t cur_start = get_start(dinucleotide_coords[current_id]);
        size_t cur_finish = get_finish(dinucleotide_coords[current_id]);
        size_t next_id = current_id + 1;
        while ((next_id) < total_d && get_start(dinucleotide_coords[next_id]) <= cur_finish + MIN_DISTANCE_COMPLEX_REGIONS) {
            cur_finish = std::max(cur_finish, get_finish(dinucleotide_coords[next_id]));
            next_id ++;

        }
        if (cur_finish > len)
            cur_finish = len;
        complex_regions.emplace_back(make_pair(cur_start, cur_finish - cur_start));
        current_id = next_id;
    }
    return complex_regions;
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
    static const size_t MAX_CONSENSUS_COVERAGE = 25;

    //Complex_regions: map from start to length;

    vector<pair<size_t, size_t>> complex_regions;
    map<size_t, vector<string>> complex_strings;
//TODO more efficient data structure

/*    vector <dinucleotide> dinucleotide_coords;
    vector <int> num_dinucleotide;
    vector<vector<size_t>> dinucleotide_read_voting;
*/

    void FillComplex() {
        complex_regions = GetComplexRegions(sequence);
        for (auto i = 0; i < complex_regions.size(); i++)
            complex_strings[complex_regions[i].first] = vector<string>();
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


    void AddRead() {

    }

    string MSAConsensus(vector<string> &s) {

//Magic consts from spoa default settings
        auto alignment_engine = spoa::AlignmentEngine::Create(
// -8 in default for third parameter(gap) opening, -6 for forth(gap extension)
                spoa::AlignmentType::kSW, 10, -8, -8, -1);  // linear gaps                                                                                                                                                          +

//TODO: process length signinficantly different with median
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
        vector<uint32_t > coverages;
        string consensus = graph.GenerateConsensus(&coverages);
        size_t pref_remove = 0;
        int suf_remove = coverages.size() - 1;
//        cout << endl;
        while (pref_remove < coverages.size() && coverages[pref_remove] < cov / 2 )
            pref_remove ++;
        while (suf_remove >= 0 && coverages[suf_remove] < cov / 2 )
            suf_remove --;
        if (pref_remove > suf_remove) {
            cout << "MSA cleanup stupid case" << endl;
            return "";
        }
        return consensus.substr(pref_remove, suf_remove - pref_remove + 1);

    }

    string checkMSAConsensus(string s, vector<string> &all) {
        size_t cons_len = s.length();
//consensus in last
        vector<size_t> all_len(all.size() - 1);
        for (size_t i = 0; i < all.size() - 1; i ++) {
            all_len[i] = all[i].length();
        }
        sort(all_len.begin(), all_len.end());
        s.erase(std::unique(s.begin(), s.end()), s.end());
        if (s.length() <= 2 * COMPLEX_EPS + 2) return "TOO SHORT";
        if (cons_len != all_len[all_len.size()/2]) {
            return "CONSENSUS LENGTH DIFFER FROM MEDIAN";
        }
        return "";
//        return "AAAAAAA!";
//What are other suspicious cases? Since we can glue two dimeric regions, commented case is  actually OK
/*        size_t len = s.length();
        for (size_t i = COMPLEX_EPS + 2; i < len - COMPLEX_EPS; i++) {
            if (s[i] != s[COMPLEX_EPS] && s[i] != s[COMPLEX_EPS + 1])
                return "On position " + to_string(i) + " of " + to_string(len) + " not dimeric nucleo";
        } */


    }


    string GenerateConsensus(Logger & logger){
        io:stringstream ss;
        vector<int> quantities(256);

        ofstream debug;
        if (debug_f != "none") {
            debug.open(debug_f + name);
            debug << name << endl;
        }
        size_t total_count = 0 ;
        size_t cur_complex_ind = 0;
#pragma omp parallel for
        for (size_t i = 0; i < complex_regions.size(); i++) {
            size_t start_pos = complex_regions[i].first;
            auto consensus = MSAConsensus(complex_strings[start_pos]);
            complex_strings[start_pos].push_back(consensus);
        }
        logger.info() << " Consenus for contig " << name << " calculated "<< endl;
        string consensus = "";

        for (size_t i = 0; i < len; ) {
            if (complex_regions.size() > 0 && complex_regions[cur_complex_ind].first == i ) {
                consensus = complex_strings[i][complex_strings[i].size() - 1];
                ss << consensus;
                auto check = checkMSAConsensus(consensus, complex_strings[i]);
 //            logger.info() << "consensus of " << complex_strings[start_pos].size() << ": " << consensus.length() << endl << "At position " <<start_pos << endl;
                if (check != ""){
                    {
                        logger.info() << "Problematic consensus starting on decompressed position " << total_count <<" " << check <<" of " <<complex_strings[i].size() - 1 << " sequences "<< endl;
                        stringstream debug_l;
                        debug_l << "lengths: ";
                        for (size_t j = 0; j < complex_strings[i].size() - 1; j++) {
                            debug_l << complex_strings[i][j].length() << " ";
                        }
                        debug_l <<" : " << consensus.length() << endl;
                        logger.info() << debug_l.str();
                         for (size_t j = 0; j < complex_strings[i].size() - 1; j++) {
                            logger.info() << complex_strings[i][j] << endl;
                        }
                        logger.info() << endl;
                        logger.info() << consensus << endl;
                    }
                }

                if (debug_f != "none" )
                    for (size_t j = 0; j < consensus.length(); j++) {
                        debug << total_count + 1 << " 0 "<< complex_strings[i].size() <<" 0 " << consensus[j] << endl;
                        total_count ++;
                    }

                i += complex_regions[cur_complex_ind].second;
                if (cur_complex_ind + 1< complex_regions.size())
                    cur_complex_ind ++;
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
//                        logger.info() << "Median disagree with average at position " << i<<" med/avg: "<< real_cov << "/" << cov << endl;
                        stringstream as;
                        for (auto j = 0; j < amounts[i][0]; j++)
                            as << int(amounts[i][j + 1]) << " ";
//                        logger.info() << as.str() << endl;
                        cov = real_cov;
                    }
                }

                quantities[quantity[i]]++;
                for (size_t j = 0; j < cov; j++) {
                    ss << sequence[i];
                    total_count++;
                }
                if (debug_f != "none") {
                    for (size_t j = 0; j < cov; j++) {
                        debug << total_count - cov + j + 1 << " " << sum[i] << " " << int(quantity[i]) << " " << cov << " "
                              << sequence[i];
                        debug << "    ";
                        for (auto jj = 0; jj < amounts[i][0]; jj++)
                            debug << int(amounts[i][jj + 1]) << " ";
                        debug << endl;
                    }
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

//We do not believe matches on the ends of match region
//TODO DO WE
    static const size_t MATCH_EPS = 1;

    static const size_t BATCH_SIZE = 10000;

    AssemblyInfo (CLParser& parser):parser(parser), logger(true, false), used_pairs(0), filtered_pairs(0), get_uncovered_only(false) {
//TODO paths;
        logging::LoggerStorage ls(".", "homopolish");
        logger.addLogFile(ls.newLoggerFile());
        logger.info() <<"Reading contigs..." << endl;
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
    bool verifyCigar(vector<cigar_pair> &cigars, int bandwidth ) {
//TODO constant??
        if (cigars.size() > 200) {
            return false;
        }
        int shift = 0;
        for(size_t i = 0; i + 1 < cigars.size(); i++) {
            if (cigars[i].type != 'M' && cigars[i+1].type != 'M') {
                return false;
            }
            if (cigars[i].type == 'D') {
                shift -= cigars[i].length;
            } else if (cigars[i].type == 'I') {
                shift += cigars[i].length;
            }
            if (shift > bandwidth || shift < -bandwidth)
                return false;
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

    std::vector<cigar_pair> getFastAln(AlignmentInfo& aln, const char * contig, const char *read) {

        size_t cur_bandwidth = SW_BANDWIDTH;
//strings, match, mismatch, gap_open, gap_extend, width
        auto cigars = align_ksw(contig, read, 1, -5, 5, 2, cur_bandwidth);
        auto str_cigars = str(cigars);
        size_t matched_l = matchedLength(cigars);
        bool valid_cigar = true;
//TODO: consts
        while ((matched_l + 300 < aln.length() || !(valid_cigar = verifyCigar(cigars, cur_bandwidth)))) {
            if (matched_l < 50) {
                logger.trace() << aln.read_id << " ultrashort alignmnent, doing nothing" << endl;
                break;
            } else {
                cur_bandwidth *= 2;
                if (cur_bandwidth > 200) {
                    break;
                }
                logger.trace() << aln.read_id << endl << str(cigars) << endl << "aln length " << aln.length()
                               << " read length " << strlen(read)
                               << " matched length " << matched_l << endl;
                cigars = align_ksw(contig, read, 1, -5, 5, 2, cur_bandwidth);
                size_t new_matched_len = matchedLength(cigars);
                logger.trace() << aln.read_id << " alignment replaced using bandwindth " << cur_bandwidth << endl
                               << str(cigars) << endl;
                if (new_matched_len <= matched_l && valid_cigar) {
                    logger.trace() << aln.read_id << " alignmnent length did not improve after moving to bandwidth " << cur_bandwidth << endl;
                    matched_l = new_matched_len;
                    break;
                }
                matched_l = new_matched_len;
            }
        }
        return cigars;
    }

    struct SplittedSequence {
        vector<string> normal_fragments;
        vector<string> complex_fragments;
        vector<size_t> normal_shifts;
        vector<size_t> complex_shifts;

        SplittedSequence(vector<pair<size_t,size_t>> &regions, string& seq) {
            size_t cur_start = 0;
            for (auto i = 0; i < regions.size(); i++) {
                normal_fragments.emplace_back(seq.substr(cur_start, regions[i].first - cur_start));
                normal_shifts.push_back(cur_start);
                complex_shifts.push_back(regions[i].first);
                cur_start = regions[i].first + regions[i].second;
                complex_fragments.emplace_back(seq.substr(regions[i].first, regions[i].second));
            }
            normal_fragments.emplace_back(seq.substr(cur_start));
            normal_shifts.push_back(cur_start);
        }
    };

    void processFragmentAlignment(AlignmentInfo& aln, const char * contig, const char *read, size_t cont_coords, size_t read_coords, vector<size_t> &quantities, ContigInfo& current_contig) {

        auto cigars = getFastAln(aln, contig, read);
        size_t matches = 0;
        size_t mismatches = 0;
        size_t indels = 0;
        size_t complex_start = -1;
        size_t contig_finish = -1;
        size_t complex_id = -1;
        size_t local_read_coords = 0;
        size_t match_regions = 0;
        size_t current_match = 0;
        for (auto it = cigars.begin(); it != cigars.end(); ++it) {
            if ((*it).type == 'M')
                match_regions ++;
        }
        for (auto it = cigars.begin(); it != cigars.end(); ++it) {
            if ((*it).type == 'I') {
                local_read_coords += (*it).length;
                indels += (*it).length;
            } else if ((*it).type == 'D') {
                cont_coords += (*it).length;
//TODO do not forget add complex regions here
                indels += (*it).length;
            } else {
                current_match ++;
                size_t cur_complex_coord = -1;
                auto complex_regions_iter = lower_bound(current_contig.complex_regions.begin(), current_contig.complex_regions.end(),make_pair(cont_coords + aln.alignment_start, size_t(0)));
                if (complex_regions_iter !=  current_contig.complex_regions.end()) {
                    cur_complex_coord = complex_regions_iter->first;
                }
//MATCH_EPS to avoid side effects near indels, zero on the ends to avoid systematic uncovered near ends of the complex regions
                size_t start_shift = MATCH_EPS;
                if (current_match == 1) start_shift = 0;
                size_t end_shift = MATCH_EPS;
                if (current_match == 1) end_shift = 0;

                for (size_t i = start_shift; i + end_shift< (*it).length; i++) {
                    size_t coord = cont_coords + aln.alignment_start + i;


//                    logger.trace() <<i <<  current_contig.sequence[coord]<< nucl(read[local_read_coords + i]) << endl;
                    if (current_contig.sequence[coord] == nucl(read[local_read_coords + i]) && current_contig.quantity[coord] != 255 && current_contig.sum[coord] < 60000) {
                        {
                            current_contig.quantity[coord]++;
                            current_contig.sum[coord] += quantities[read_coords + local_read_coords+ i];
                            if (current_contig.amounts[coord][0] < ContigInfo::VOTES_STORED)
#pragma omp critical
                                current_contig.amounts[coord][++current_contig.amounts[coord][0]] = quantities[read_coords + local_read_coords + i];
                        }
                        matches ++;
                    } else {
                        mismatches ++;
                    }
//Only complete traversion of complex regions taken in account;
//                    if (current_contig.complex_regions.find(coord) != current_contig.complex_regions.end()) {
                    if (coord == cur_complex_coord) {
//                        logger.trace() << "entered complex at coord" << cur_complex_coord;
                        size_t complex_len = complex_regions_iter->second;
                        if (local_read_coords + complex_len < matchedLength(cigars)) {
                            complex_start = local_read_coords + i;
                            contig_finish = coord + complex_len;
                            complex_id = coord;
                        }
                        complex_regions_iter ++;
                        if (complex_regions_iter != current_contig.complex_regions.end())
                            cur_complex_coord = complex_regions_iter->first;
                    }
                    if (coord == contig_finish) {
                        contig_finish = -1;
                    }
                }
                local_read_coords += (*it).length;
                cont_coords += (*it).length;
            }
        }
        if (matches < mismatches * 3)
            logger.info()<<"For read too much mismatches " << aln.read_id << " matches/MM: " << matches << "/" << mismatches << endl;
    }

    void processReadPair (string& read, AlignmentInfo& aln) {
//        logger.info() << read.id << endl;
        Sequence read_seq (read);
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

        std:vector<size_t> quantities;
        quantities.resize(compressed_read.length());

        int cur_ind = 0;
        for (size_t i = 0; i < compressed_read.size(); i++){
            size_t count = 0;
            while (cur_ind < read_seq.size() && nucl(read_seq[cur_ind]) == compressed_read[i]) {
                count ++;
                cur_ind ++;
            }
            quantities[i] = count;
        }

        auto read_complex_regions = GetComplexRegions(compressed_read);
        vector<pair<size_t, size_t>> contig_complex_regions = GetComplexRegions(contig_seq);
//Do not want to take them from global contigs' complex regions because of side effects
/*
        auto complex_regions_iter = lower_bound(current_contig.complex_regions.begin(), current_contig.complex_regions.end(),make_pair(aln.alignment_start, size_t(0)));
        while(complex_regions_iter !=  current_contig.complex_regions.end()) {
            //complex region completely inside alignment, NEED TO check start also
            if (complex_regions_iter->first + complex_regions_iter->second < aln.alignment_end) {
                contig_complex_regions.push_back(*complex_regions_iter);
                complex_regions_iter ++;
            } else {
                break;
            }
        }
        */
        size_t contig_complex_size = contig_complex_regions.size();
        size_t read_complex_size = read_complex_regions.size();

        if (contig_complex_regions.size() != read_complex_regions.size()) {
            logger.info() << aln.read_id << " different complex regions amounts " << contig_complex_size << " vs " << read_complex_size << endl;
            return;
        }


        SplittedSequence read_splt(read_complex_regions, compressed_read);
        SplittedSequence contig_splt(contig_complex_regions, contig_seq);
        for (size_t i = 0; i < read_splt.normal_fragments.size(); i ++) {
            logger.trace() << aln.read_id << "frag  " << i << " start " << contig_splt.normal_shifts[i] + aln.alignment_start << endl;
//TODO  HERE SHOULD BE SIZE CHECK
            {
//                logger.trace() << contig_splt.normal_fragments[i] << endl << read_splt.normal_fragments[i] << endl;
                processFragmentAlignment(aln, contig_splt.normal_fragments[i].c_str(),
                                         read_splt.normal_fragments[i].c_str(), contig_splt.normal_shifts[i],
                                         read_splt.normal_shifts[i], quantities, current_contig);
            }
        }
//        logger.info() << "normal fragments processed";
        for (auto i = 0; i < read_splt.complex_fragments.size(); i++) {
//TODO: check on sizes;
            if (contig_splt.complex_shifts[i] < MIN_DISTANCE_COMPLEX_REGIONS ||  contig_splt.complex_shifts[i] +  contig_splt.complex_fragments[i].length() + MIN_DISTANCE_COMPLEX_REGIONS > contig_seq.length()) {
                logger.info() << "For read " <<aln.read_id <<" skipping region " << i <<" (of "<<contig_splt.complex_fragments.size() << " end of contig's fragment" << endl;
                continue;
            }
            if (read_splt.complex_shifts[i] < MIN_DISTANCE_COMPLEX_REGIONS ||  read_splt.complex_shifts[i] +  read_splt.complex_fragments[i].length() + MIN_DISTANCE_COMPLEX_REGIONS > compressed_read.length()) {
                logger.info() << "For read " <<aln.read_id <<" skipping region " << i <<" (of "<<read_splt.complex_fragments.size() << " end of read" << endl;
                continue;
            }
            size_t pos = contig_splt.complex_shifts[i] + aln.alignment_start;
            if (current_contig.complex_strings.find(pos) != current_contig.complex_strings.end())
                current_contig.complex_strings[pos].push_back(uncompress(read_splt.complex_shifts[i], read_splt.complex_shifts[i] + read_splt.complex_fragments[i].length(), read_seq.str()));
            else {
                logger.info() << "No complex region at position " << pos << endl;
            }
        }
//        logger.info() << "complex fragments processed" << endl;

//        logger.info() << "Matches "<< matches << " mismatches " << mismatches <<" indels " << indels << endl;
        return;
    }
    void removeWhitespace(string & s){
        auto white = s.find(' ');
        if (white != string::npos) {
            s = s.substr(0, white);
        }
    }

    void processBatch(vector<string>& contigs, vector<AlignmentInfo>& alignments){
        size_t len = contigs.size();
#pragma omp parallel for
        for (size_t i = 0; i < len; i++) {
            processReadPair(contigs[i], alignments[i]);
        }
    }

    void process() {

        ifstream compressed_reads;
        ofstream corrected_contigs;
        corrected_contigs.open(parser.getValue("output"));
        compressed_reads.open(parser.getValue("alignments"));
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
        vector<AlignmentInfo> align_batch;
        vector<string> contig_batch;
        while (!compressed_reads.eof()) {
            bool reads_over = false;
            while (cur.id != cur_compressed) {

                if (reader.eof()) {
                    logger.info() << "Reads are over\n";
                    reads_over = true;
                    break;
                }
                cur = reader.read();
//Some tools clip read names after whitespaces
                removeWhitespace(cur.id);
                reads_count ++;
                if (reads_count % 1000 == 0) {
                    logger.info() << "Processed " << reads_count << " original reads " << endl;
                }
            }
            if (reads_over) {
                break;
            }
            align_batch.push_back(cur_align);
            contig_batch.push_back(cur.seq);
//            processReadPair(cur, cur_align);
//TODO:: appropriate logic for multiple alignment
            do {
                cur_align = readAlignment(compressed_reads);
                aln_count ++;
                if (aln_count % BATCH_SIZE == 0) {
                    logger.info() << "Batch of size " <<BATCH_SIZE <<" created, processing" << endl;
                    processBatch(contig_batch, align_batch);

                    logger.info() << "Processed " << aln_count << " compressed mappings " << endl;
                    contig_batch.resize(0);
                    align_batch.resize(0);
                    //exit(0);
                }

                if (cur_compressed == cur_align.read_id) {
                    align_batch.push_back(cur_align);
                    contig_batch.push_back(cur.seq);

//                    contig_batch.emplace_back(cur);
//                    processReadPair(cur, cur_align);
                }
            } while (cur_compressed == cur_align.read_id);
/*            if (aln_count %2000 == 1999)
                break; */
            cur_compressed = cur_align.read_id;

        }
        processBatch(contig_batch, align_batch);
        logger.info() << "Processed final batch of " << align_batch.size() << " compressed reads " << endl;

        logger.info() << "Reads processed, used " <<used_pairs << " filtered by compressed length " << filtered_pairs << endl;
        for (auto& contig: contigs){
            logger.info() << "Generating consensus for contig " << contig.first << endl;
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
    CLParser parser({"alignments=", "contigs=", "output=", "debug=none", "threads="}, {"reads"},
                    {});

    parser.parseCL(argc, argv);

    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters" << std::endl;
        std::cout << parser.check() << std::endl;
        return 1;
    }

    omp_set_num_threads(stoi(parser.getValue("threads")));
    AssemblyInfo assemblyInfo(parser);
    assemblyInfo.process();

 }
