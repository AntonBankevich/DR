#include <common/verify.hpp>
#include <sequences/seqio.hpp>
#include <common/cl_parser.hpp>

int main(int argc, char **argv) {
    CLParser parser({"dimer-compress=1000000000,1000000000,1", "contigs=", "bed="}, {}, {},"");
    parser.parseCL(argc, argv);
    StringContig::homopolymer_compressing = true;
    StringContig::SetDimerParameters(parser.getValue("dimer-compress"));

    std::ifstream is;
    is.open(parser.getValue("bed"));
    io::Library contigs_lib = {std::experimental::filesystem::path(parser.getValue("contigs"))};
    io::SeqReader reader(contigs_lib);
    std::unordered_map<std::string, std::string> contigs;
    for(StringContig s : reader) {
        contigs.emplace(s.id, std::move(s.seq));
    }
    std::string line;
    while(getline(is,line)) {
        std::vector<std::string> s = split(line);
        if(s.size() != 3)
            break;
        size_t start = std::stoull(s[1]);
        size_t end = std::stoull(s[2]);
        std::string &seq = contigs[s[0]];
        start = StringContig("c1", seq.substr(0, start)).makeContig().size();
        end = StringContig("c2", seq.substr(0, end)).makeContig().size();
        std::cout << s[0] << "\t" << start << "\t" << end << "\n";
    }
    is.close();
}
