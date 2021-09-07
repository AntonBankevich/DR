#include <common/verify.hpp>
#include <sequences/seqio.hpp>
#include <common/cl_parser.hpp>

class RawSegment {
public:
    std::string id = "";
    size_t left = 0;
    size_t right = 0;
    RawSegment() = default;
    RawSegment(const std::string &_id, size_t _left, size_t _right) : id(_id), left(_left), right(_right){}

    bool valid() const {
        return !id.empty();
    }

    RawSegment cap(const RawSegment &other) const {
        if (id == other.id && std::max(left, other.left) < std::min(right, other.right)) {
            return {id, std::max(left, other.left), std::min(right, other.right)};
        } else {
            return {};
        }
    }

    size_t size() const {
        return right - left;
    }
};

std::vector<RawSegment> collectSegments(const std::vector<std::string> &files) {
    std::vector<RawSegment> interest;
    for(const std::string &p : files) {
        std::ifstream is;
        is.open(p);
        std::string line;
        while(getline(is,line)) {
            std::vector<std::string> s = split(line);
            if(s.size() != 3)
                break;
            size_t start = std::stoull(s[1]);
            size_t end = std::stoull(s[2]);
            interest.emplace_back(s[0], start, end);
        }
        is.close();
    }
    return interest;
}

bool compare(const RawSegment &seg1, const RawSegment &seg2) {
    if(seg1.size() != seg2.size())
        return seg1.size() > seg2.size();
    if(seg1.id != seg2.id)
        return seg1.id < seg2.id;
    if(seg1.left != seg2.left)
        return seg1.left < seg2.left;
    return seg1.right < seg2.right;
}

int main(int argc, char **argv) {
    CLParser parser({"quast="}, {"bed"}, {},"");
    parser.parseCL(argc, argv);
    std::vector<RawSegment> interest = collectSegments(parser.getListValue("bed"));
    std::experimental::filesystem::path dir(parser.getValue("quast"));
    for(const auto &p : std::experimental::filesystem::directory_iterator(dir / "contigs_reports")){
        std::string fname = p.path().filename().string();
        std::string prefix = "all_alignments_";
        std::string suffix = ".tsv";
        std::vector<std::vector<RawSegment>> segs(interest.size());
        if (startsWith(fname, prefix) && endsWith(fname, suffix)) {
            std::string name = fname.substr(prefix.size(), fname.size() - prefix.size() - suffix.size());
            std::string line;
            std::ifstream is;
            is.open(p.path());
            getline(is,line);
            size_t left = size_t(-1);
            size_t right = 0;
            std::string chr = "";
            while(true) {
                if(!getline(is,line)) {
                    break;
                }
                std::vector<std::string> s = split(line);
                left = std::min<size_t>(left, std::stoull(s[0]));
                right = std::max<size_t>(right, std::stoull(s[1]));
                VERIFY(chr.empty() || chr == s[4]);
                chr = s[4];
                getline(is,line);
                if(startsWith(line, "CONTIG") || startsWith(line, "relocation,") ||
                            startsWith(line, "translocation") || startsWith(line, "invert")) {
                    RawSegment newSeg(chr, left, right);
                    for(size_t i = 0; i < interest.size(); i++) {
                        if(interest[i].cap(newSeg).valid())
                            segs[i].emplace_back(interest[i].cap(newSeg));
                    }
                    left = size_t(-1);
                    right = 0;
                    chr = "";
                }
            }
            is.close();
            for(size_t i = 0; i < interest.size(); i++) {
                std::sort(segs[i].begin(), segs[i].end(), compare);
                VERIFY(segs[i].size() == 0 || segs[i].front().size() >= segs[i].back().size());
                size_t len = interest[i].size();
                size_t tmp = 0;
                size_t j = 0;
                for(; j < segs[i].size() && tmp < len * 95 / 100; j++)
                    tmp += segs[i][j].size();
                if(tmp >= len * 95 / 100) {
                    std::cout << name << "\t" << interest[i].id << "\t" << interest[i].left << "\t" << interest[i].right << "\t" << j << std::endl;
                } else {
                    std::cout << name << "\t" << interest[i].id << "\t" << interest[i].left << "\t" << interest[i].right << "\t" << "NA" << std::endl;
                }
            }
        }
    }
}