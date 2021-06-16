#include <sequences/verify.hpp>
#include <sequences/seqio.hpp>

int main(int argc, char **argv) {
    VERIFY(argc == 2);
    StringContig::needs_compressing = false;
    io::Library libReads = {std::experimental::filesystem::path(argv[1])};
    io::SeqReader reader(libReads);
    std::vector<std::pair<Sequence, Contig>> contigs;
    for(StringContig s : reader) {
        Contig c = s.makeContig();
        if(c.seq <= !c.seq)
            std::cout << ">" << c.id << "\n" << c.seq << "\n";
    }
}