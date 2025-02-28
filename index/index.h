#ifndef INDEX_H
#define INDEX_H

#include <bits/stdc++.h>
#include <chrono>
#include <algorithm>

#include <sdsl/suffix_trees.hpp>
#include <sdsl/rmq_support.hpp>


struct FastaData {
    std::vector<std::string> names;
    std::vector<std::string> sequences;
};

struct Indexcfg{
    bool rebuild = false;
    std::filesystem::path text_file;
    std::filesystem::path base_path; 
};

struct mem_occ{
    int index;
    int length;
    int number;
    std::vector<int> locations;
    std::vector<int> introns;
    mem_occ(int index,int length, int number, std::vector<int> locations,std::vector<int> introns) :index(index), length(length), number(number),locations(locations), introns(introns){}
};

class Index
{
private:
    // typedef sdsl::cst_sada<sdsl::csa_wt<>, sdsl::lcp_support_sada<>> cst_type;
    typedef sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector_il<>>, 1024, 1<<20>, sdsl::lcp_support_sada<>> cst_type;
    bool rebuild = false;
    sdsl::csa_bitcompressed<> tmp_csa;

public:
    typedef cst_type::size_type size_type;

    Indexcfg config;
    size_t text_size;
    sdsl::rank_support_v<> rankB;
    sdsl::bit_vector B;
    cst_type cst;

    sdsl::rmq_succinct_sct<true> rmq_sa_min;
    sdsl::rmq_succinct_sct<false> rmq_sa_max;

    double index_size;
    std::vector<mem_occ> occurences;

    Index(std::filesystem::path output_path);
    ~Index();

    double build(); //  construct fm index supporting mem searching and rmq structure over SA
    double locate(std::string pattern);            //  locate all MeMs
};
#endif //  INDEX_H
