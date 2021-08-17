#include <common/verify.hpp>
#include <experimental/filesystem>
#include <fstream>
#include <unistd.h>
#include <common/dir_utils.hpp>
#include <common/string_utils.hpp>

std::string getFirst(const std::experimental::filesystem::path &cfile) {
    std::string res;
    std::ifstream is;
    is.open(cfile);
    std::getline(is, res);
    is.close();
    return std::move(res);
}

void removeFirst(const std::experimental::filesystem::path &cfile) {
    std::vector<std::string> contents;
    std::string next;
    std::ifstream is;
    is.open(cfile);
    while(std::getline(is, next)) {
        if(!next.empty()) {
            contents.emplace_back(next);
        }
    }
    is.close();
    std::ofstream os;
    os.open(cfile);
    contents.erase(contents.begin());
    for(std::string &line : contents) {
        os << line << "\n";
    }
    os.close();
}
void append(const std::experimental::filesystem::path &cfile, const std::string &next) {
    std::ofstream os;
    os.open(cfile, std::ios_base::app);
    os << next << "\n";
    os.close();
}

int main(int argc, char **argv) {
    VERIFY(argc == 3);
    std::experimental::filesystem::path cfile(argv[1]);
    std::experimental::filesystem::path dir(argv[2]);
    ensure_dir_existance(dir);
    ensure_dir_existance(dir / "logs");
    std::experimental::filesystem::path log_path = dir/"log.txt";
    size_t max = 0;
    for(const std::experimental::filesystem::path &p :std::experimental::filesystem::directory_iterator(dir / "logs")){
        max = std::max<size_t>(max, std::stoull(p.filename()));
    }
    bool slept = false;
#pragma clang diagnostic push
#pragma ide diagnostic ignored "EndlessLoop"
    while(true) {
        std::string command = getFirst(cfile);
        if(command.empty()) {
            if(!slept) {
                std::cout << "No command. Sleeping." << std::endl;
                slept = true;
            }
            sleep(60);
        } else {
            slept = false;
            std::cout << "Running new command: " << command << std::endl;
            std::experimental::filesystem::path log(dir/"logs"/itos(max + 1));
            std::cout << "Printing output to file " << log << std::endl;
            std::ofstream os;
            os.open(log_path, std::ios_base::app);
            os << (max + 1) << " started\n" << command << "\n";
            os.close();
            int res = system((command + " > " + log.string() + " 2>&1").c_str());
            std::cout << "Finished running command: " << command << std::endl;
            std::cout << "Return code: " << res << "\n\n" << std::endl;
            std::ofstream os1;
            os1.open(log_path, std::ios_base::app);
            os1 << (max + 1) << " finished " << res << "\n\n";
            os1.close();
            max++;
            removeFirst(cfile);
            append(dir / "finished.txt", command);
        }
    }
#pragma clang diagnostic pop
}