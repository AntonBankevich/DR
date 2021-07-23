#include <common/verify.hpp>
#include <sequences/seqio.hpp>

int main(int argc, char **argv) {
    VERIFY(argc == 2);
    StringContig::homopolymer_compressing = false;
    std::experimental::filesystem::path dir(argv[1]);
    for(const auto &f : std::experimental::filesystem::directory_iterator(dir)) {
        if(!endsWith(f.path().filename().string(), ".fasta"))
            continue;
        std::string name = f.path().filename().string().substr(0, f.path().filename().string().size() - 6);
        io::SeqReader reader({f.path()});
        for(StringContig s : reader) {
            Contig c = s.makeContig();
            if(c.seq <= !c.seq)
                std::cout << ">" << name << "." << c.id << "\n" << c.seq << "\n";
        }
    }
}