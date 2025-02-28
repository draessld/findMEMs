#include <filesystem>
#include <string>
#include <iostream>

#include "index/index.h"

std::string usage = "./locate <input_file> <pattern>";

using namespace std;

template <typename T>
void print_MEMs(std::vector<T> occurences, std::string pattern)
{
    for (size_t i = 0; i < occurences.size(); i++)
    {
        std::cout << pattern.substr(occurences[i].index,occurences[i].length)<<":"<< occurences[i].number <<"\t";
        for (size_t j = 0; j < occurences[i].locations.size(); j++)
        {
            std::cout << occurences[i].locations[j] << "{" << occurences[i].introns[j]<<"}\t";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

int main(int argc, char **argv)
{

    /*  Handle input parameters */
    if (argc != 3)
    {
        cerr << "Bad arguments\n " << usage;
        return -1;
    }
    
    std::filesystem::path input_path = argv[1];
    std::string pattern_str = argv[2];
    std::vector<std::string> patterns;
    std::string tmp;
    

    std::ifstream ifs(pattern_str);
    if(ifs.is_open()){
        // Reading patten file
        tmp = std::filesystem::path(pattern_str).extension();
        if(tmp == ".fa" || tmp == ".fasta")
        {
            //  every fasta record is pattern
            while (std::getline(ifs, tmp))
            {
                if (tmp.empty() || tmp[0] == '>')
                    continue;
                patterns.push_back(tmp.substr(0, tmp.size() - 1));
            }
        }else{
            //  every line one pattern
            while (std::getline(ifs, tmp))
            {
                if (tmp.empty())
                    continue;
                patterns.push_back(tmp.substr(0, tmp.size() - 1));
            }
        }
    }else{
        // input as a string
        patterns.push_back(pattern_str);
    }

    //  Load index
    Index index = Index(input_path);
    index.build();

    //  locate pattern
    for (auto pattern : patterns)
    {
        std::cout << '>' << pattern << '\t' << index.locate(pattern) << "us" << '\t' << index.occurences.size() << std::endl;
        print_MEMs(index.occurences, pattern);
        index.occurences.clear();
    }

    return 0;
}