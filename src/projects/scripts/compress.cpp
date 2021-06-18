#include <sequences/verify.hpp>
#include <sequences/seqio.hpp>

int main(int argc, char **argv) {
    VERIFY(argc == 2);
    StringContig::homopolymer_compressing = true;
    io::Library libReads = {std::experimental::filesystem::path(argv[1])};
    io::SeqReader reader(libReads);
    for(StringContig s : reader) {
        Contig c = s.makeContig();
        std::cout << ">" << c.id << "\n" << c.seq << "\n";
    }
}
