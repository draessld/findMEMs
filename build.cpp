#include <filesystem>
#include <string>
#include <iostream>

#include "index/index.h"

using namespace std;

std::string usage = "./build <input_file>";

int main(int argc, char **argv)
{
    /*  Handle input parameters */
    if (argc != 2)
    {
        cerr << "Bad arguments\n " << usage;
        return -1;
    }
    
    std::filesystem::path input_path = argv[1];

    //  create index
    Index index = Index(input_path); //  load or build index
    double build_result = index.build();

    std::cout << "Build time: " << build_result << "μs" << std::endl;
    return 0;
}