#include <common/verify.hpp>
#include <sequences/seqio.hpp>

int main(int argc, char **argv) {
    VERIFY(argc == 2);
    StringContig::homopolymer_compressing = false;
    io::Library libReads = {std::experimental::filesystem::path(argv[1])};
    io::SeqReader reader(libReads);
    std::unordered_map<std::string, Sequence> m;
    for(StringContig s : reader) {
        Contig c = s.makeContig();
        m[c.getId()] = c.seq;
    }
    for(auto &it : m) {
        std::cout << ">" << it.first << "\n" << it.second << "\n";
    }
}