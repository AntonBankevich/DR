#include <sequences/verify.hpp>
#include <sequences/seqio.hpp>

int main(int argc, char **argv) {
    VERIFY(argc == 1);
    StringContig::needs_compressing = true;
    io::Library libReads = {std::experimental::filesystem::path(argv[0])};
    io::SeqReader reader(libReads);
    for(StringContig s : reader) {
        Contig c = s.makeContig();
        std::cout << c.id << "\n" << c.seq << "\n";
    }
}
