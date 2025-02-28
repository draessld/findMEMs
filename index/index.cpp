#include "index.h"

Index::Index(std::filesystem::path text_file)
{
    std::string filename = text_file.filename();
    std::string ext = text_file.extension().c_str();
    config.base_path = text_file.replace_extension(ext + ".index");
    if (!std::filesystem::exists(config.base_path))
        std::filesystem::create_directories(config.base_path);
    config.base_path.append(filename);
    config.text_file = text_file.replace_filename(filename);

    if (ext == ".fa" || ext == ".fasta")
    {
        //  Read fasta
        std::cout << "Reading fasta file" << std::endl;
        std::ifstream ifs(config.text_file);

        config.text_file = config.base_path.replace_extension(".txt");
        std::ofstream ofs(config.text_file);
        std::string line;
        if (!ifs.is_open())
        {
            std::cerr << "Error opening file: " << filename << std::endl;
        }

        std::getline(ifs, line);
        while (std::getline(ifs, line))
        {
            if (line.empty())
                continue;
            if (line[0] == '>')
                ofs << '#';
            else
                ofs << line.substr(0, line.size() - 1);
        }

        ifs.close();
        ofs.close();
        config.base_path.replace_extension("");
    }

    std::cout << "Index destination on " << config.base_path << std::endl;
    std::cout << "Text destination on " << config.text_file << std::endl;
}

Index::~Index() {}

double Index::build()
{

    std::cout << "-=-=-=-=-=-   Building index   ...   " << std::endl;

    auto base = std::chrono::high_resolution_clock::now();
    auto startTime = base;

    //  build compressed suffix tree
    if (!sdsl::load_from_file(cst, config.base_path.replace_extension(".cst")))
    {
        std::ifstream in(config.text_file);
        if (!in)
        {
            std::cout << "ERROR: File " << config.text_file << " does not exist. Exit." << std::endl;
        }
        std::cout << "No compact suffix tree on " << config.base_path.replace_extension(".cst") << " located. Building now ...   ";
        sdsl::construct(cst, config.text_file, 1);
        sdsl::store_to_file(cst, config.base_path.replace_extension(".cst"));
        std::cout << " ==> DONE in " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startTime).count() << "um" << std::endl;
    }
    startTime = std::chrono::high_resolution_clock::now();

    // //  build RMaxQ SA support
    // if (!sdsl::load_from_file(rmq_sa_max, config.base_path.replace_extension(".rmq_max")))
    // {
    //     std::cout << "No RMQ for SA on " << config.base_path.replace_extension(".rmq_max") << " located. Building now ... ";

    //     sdsl::construct(tmp_csa, config.text_file, 1);
    //     sdsl::util::assign(rmq_sa_max, sdsl::rmq_succinct_sct<false>(&tmp_csa));
    //     sdsl::store_to_file(rmq_sa_max, config.base_path.replace_extension(".rmq_max"));

    //     std::cout << " ==> DONE in " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startTime).count() << "um" << std::endl;
    // }
    // startTime = std::chrono::high_resolution_clock::now();

    // //  build RMinQ SA support
    // if (!sdsl::load_from_file(rmq_sa_min, config.base_path.replace_extension(".rmq_min")))
    // {
    //     std::cout << "No RMQ for SA on " << config.base_path.replace_extension(".rmq_min") << " located. Building now ... ";
    //     sdsl::construct(tmp_csa, config.text_file, 1);
    //     sdsl::util::assign(rmq_sa_min, sdsl::rmq_succinct_sct<true>(&tmp_csa));
    //     sdsl::store_to_file(rmq_sa_min, config.base_path.replace_extension(".rmq_min"));
    //     std::cout << " ==> DONE in " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startTime).count() << "um" << std::endl;
    // }
    // startTime = std::chrono::high_resolution_clock::now();

    text_size = cst.size();

    //  build rank support for $ character
    if (!sdsl::load_from_file(B, config.base_path.replace_extension("ranksupp")))
    {
        std::ifstream in(config.text_file);
        if (!in)
        {
            // std::cout << "ERROR: File " << config.text_file << " does not exist. Exit." << std::endl;
        }
        std::cout << "No rank support  on " << config.base_path.replace_extension("ranksupp") << " located. Building now ...  ";
        size_t i = 0;
        B.resize(text_size);
        char c;
        while (in.get(c))
        {
            B[i] = 0;
            if (c == '#')
                B[i++] = 1;
            else
                i++;
        }
        sdsl::store_to_file(B, config.base_path.replace_extension("ranksupp"));
        std::cout << " ==> DONE in " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startTime).count() << "um" << std::endl;
    }
    sdsl::util::assign(rankB, &B);

    double index_size = size_in_mega_bytes(cst) + size_in_mega_bytes(B) + size_in_mega_bytes(rmq_sa_min) + size_in_mega_bytes(rmq_sa_max);
    std::cout << "-=-=-=-=-=-   Building index - DONE, size:  " << index_size << " MiB." << std::endl;
    // return index_size;
    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - base).count();
}

double Index::locate(std::string pattern)
{
    auto base = std::chrono::high_resolution_clock::now();

    auto r = cst.root();
    size_t length = 1;
    auto o = r;

    for (int i = pattern.size() - 1; i >= 0; i--) //  at most 2*pattern.size() iterations
    {
        //  propagate into tree
        r = o;
        o = cst.wl(r, pattern[i]);

        if (cst.depth(o) == 0)
        {
            //  character was not found - returns root

            //  report MEM
            if (length > 1)
            {
                std::vector<int> locations;
                std::vector<int> introns;
                for (size_t x = cst.lb(r); x <= cst.rb(r); ++x)
                {
                    locations.push_back(cst.csa[x]);
                    introns.push_back(rankB(cst.csa[x]));
                }
                occurences.emplace_back(i + 1, length - 1,locations.size(),locations,introns);

            }
            if (r != o) // no change from previous = both of them must be roots
                i++;

            while (cst.depth(cst.wl(r, pattern[i - 1])) == 0)
            {
                r = cst.parent(r);
                // std::cout << "." << cst.depth(r) << (cst.depth(r)==0) << std::endl;

                if (cst.depth(r) == 0)
                    break;
            }
            length = cst.depth(r);
            o = r;
        }
        if (i == 0)
        {
            //  end of the pattern, report MEM
            if (length != 0)
            {
                
                std::vector<int> locations;
                std::vector<int> introns;
                for (size_t x = cst.lb(r); x <= cst.rb(r); ++x)
                {
                    locations.push_back(cst.csa[x]);
                    introns.push_back(rankB(cst.csa[x]));
                }
                occurences.emplace_back(i, length, locations.size(), locations, introns);
            }
            break;
        }
        length++;
    }

    // return

    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - base).count();
}