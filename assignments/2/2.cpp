#include <omp.h>  // used here just for omp_get_wtime()
#include <cstring>
#include <vector>
#include <set>
#include <string>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>

using umap=std::unordered_map<std::string, uint64_t>;
using pair=std::pair<std::string, uint64_t>;
struct Comp {
	bool operator ()(const pair& p1, const pair& p2) const {
		return p1.second > p2.second;
	}
};
using ranking=std::multiset<pair, Comp>;

uint64_t total_words{0};
volatile uint64_t extraworkXline{0};

void printMap(umap myMap) {
    printf("MAP:\n");
    for(auto& elem : myMap) {
        printf("%s, %ld\n", elem.first.c_str(), elem.second);
    }
}

void reduce_umaps(umap& output, umap& input)
{
	if(!output.empty() || !input.empty()){
		// int i = omp_get_thread_num();
		// int n = omp_get_num_threads();
		// printf("reduce_umaps %d/%d\n", i,n);
		// printMap(output);
		// printMap(input);
		for(auto& X : input) {
			// printf("first %s\n",X.first.c_str());
			// printf("second %ld\n", X.second);
			
			auto t = output.find(X.first);
			if (t == output.end()) {
				output.insert(pair{X.first, X.second});
			} else {
				output.at(X.first) += X.second; //Will throw if X.first doesn't exist in output. 
			}
		}
	}
}

#pragma omp declare reduction(umap_reduction : \
   umap : reduce_umaps(omp_out, omp_in)) \
    initializer(omp_priv(omp_orig))
    
void tokenize_line(const std::string& line, umap& UM) {
    char *tmpstr;
    char *token = strtok_r(const_cast<char*>(line.c_str()), " \r\n", &tmpstr);
	// std::cout<<"extraworkXline: "<<extraworkXline<<std::endl;
	
	#pragma omp parallel num_threads(extraworkXline) firstprivate(token) shared(tmpstr)
    {	
		// int i = omp_get_thread_num();
        // int n = omp_get_num_threads();
        while(token) {
			// if(token)
            // 	printf("hello from thread %d of %d with token %s\n", i, n, token);
        	// else
            // 	printf("toke is null\n");
            ++UM[std::string(token)];
			
			#pragma omp critical
			{
		    	++total_words;
			}
            
			token = strtok_r(NULL, " \r\n", &tmpstr);
        }
    }

    // #pragma omp parallel for schedule(dynamic) firstprivate(token) shared(tmpstr) reduction(umap_reduction:UM) reduction(+:total_words)
    // for(volatile uint64_t j = 0; j<extraworkXline; j++) {
    //     while(token) {
    //         int i = omp_get_thread_num();
    //         int n = omp_get_num_threads();
    //         if(token)
    //             printf("hello from thread %d of %d with token %s\n", i, n, token);
    //         else
    //             printf("toke is null\n");
                 
    //         if(token) {
    //             ++UM[std::string(token)];
	// 			++total_words;
	// 			token = strtok_r(NULL, " \r\n", &tmpstr);
    //         }
    //     }
    // }
	//for(volatile uint64_t j{0}; j<extraworkXline; j++);
}
std::vector<std::string> readLinesParallel(const std::string& filename) {
    std::vector<std::string> lines;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return lines;
    }

    // Get the number of lines in the file
    file.seekg(0, std::ios::end);
    size_t fileSize = file.tellg();
    file.seekg(0, std::ios::beg);
    size_t numLines = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
    file.seekg(0, std::ios::beg);

    // Reserve space for the lines
    lines.reserve(numLines);

    // Read lines in parallel
    #pragma omp parallel
    {
        // #pragma omp single
        // {
        //     std::cout << "Reading file with " << omp_get_num_threads() << " threads." << std::endl;
        // }

        while (std::getline(file, line)) {
            #pragma omp task
            {
                lines.push_back(line);
            }
        }
    }

    file.close();

    return lines;
}

std::vector<std::string> getVectorLines(const std::string& filename) {
	std::ifstream file(filename, std::ios_base::in);
	std::vector<std::string> returnV;
	if (file.is_open()) {
		std::string line;
		while(std::getline(file, line)) {
			if (!line.empty()) {
				returnV.emplace_back(line);
			}
		}
	}
	file.close();
	return returnV;
}

void compute_fileReadParallel(const std::string& filename, umap& UM) {
	std::vector<std::string> vectorLines;
	#pragma omp parallel private(vectorLines)
	{
		vectorLines = readLinesParallel(filename);
		//std::cout<<vectorLines.size()<<std::endl;
		#pragma omp parallel for schedule(dynamic) firstprivate(vectorLines)
		for (auto line : vectorLines) { // forse vanno utilizzati il numero dei thread anche qui degli extra workers
			//printf("line %s\n", line.c_str());
			tokenize_line(line, UM);
		}
	}
	
}

void compute_file(const std::string& filename, umap& UM) {
	std::vector<std::string> vectorLines = getVectorLines(filename);
	//std::cout<<vectorLines.size()<<std::endl;
	#pragma omp parallel for schedule(dynamic) firstprivate(vectorLines)
	for (auto line : vectorLines) { // forse vanno utilizzati il numero dei thread anche qui degli extra workers
		//printf("line %s\n", line.c_str());
		tokenize_line(line, UM);
	}
}

void compute_fileOld(const std::string& filename, umap& UM) {
	std::ifstream file(filename, std::ios_base::in);
	if (file.is_open()) {
		std::string line;
		std::vector<std::string> V;
        while(std::getline(file, line)) {
			#pragma omp task firstprivate(line)
			{
				if (!line.empty()) {
					//printf("line %s\n", line.c_str());
					tokenize_line(line, UM);
				}
			}
		}
	} 
	file.close();
}
int main(int argc, char *argv[]) {

	auto usage_and_exit = [argv]() {
		std::printf("use: %s filelist.txt [extraworkXline] [topk] [showresults]\n", argv[0]);
		std::printf("     filelist.txt contains one txt filename per line\n");
		std::printf("     extraworkXline is the extra work done for each line, it is an integer value whose default is 0\n");
		std::printf("     topk is an integer number, its default value is 10 (top 10 words)\n");
		std::printf("     showresults is 0 or 1, if 1 the output is shown on the standard output\n\n");
		exit(-1);
	};

	std::vector<std::string> filenames;
	size_t topk = 10;
	bool showresults=false;
	if (argc < 2 || argc > 5) {
		usage_and_exit();
	}

	if (argc > 2) {
		try { extraworkXline=std::stoul(argv[2]);
		} catch(std::invalid_argument const& ex) {
			std::printf("%s is an invalid number (%s)\n", argv[2], ex.what());
			return -1;
		}
		if (argc > 3) {
			try { topk=std::stoul(argv[3]);
			} catch(std::invalid_argument const& ex) {
				std::printf("%s is an invalid number (%s)\n", argv[3], ex.what());
				return -1;
			}
			if (topk==0) {
				std::printf("%s must be a positive integer\n", argv[3]);
				return -1;
			}
			if (argc == 5) {
				int tmp;
				try { tmp=std::stol(argv[4]);
				} catch(std::invalid_argument const& ex) {
					std::printf("%s is an invalid number (%s)\n", argv[4], ex.what());
					return -1;
				}
				if (tmp == 1) showresults = true;
			}
		}
	}
	
	if (std::filesystem::is_regular_file(argv[1])) {
		std::ifstream file(argv[1], std::ios_base::in);
		if (file.is_open()) {
			std::string line;
			while(std::getline(file, line)) {
				if (std::filesystem::is_regular_file(line))
					filenames.push_back(line);
				else
					std::cout << line << " is not a regular file, skipt it\n";
			}					
		} else {
			std::printf("ERROR: opening file %s\n", argv[1]);
			return -1;
		}
		file.close();
	} else {
		std::printf("%s is not a regular file\n", argv[1]);
		usage_and_exit();
	}

	// used for storing results
	umap UM;

	// start the time
	auto start = omp_get_wtime();	

	#pragma omp parallel for schedule(dynamic) reduction(umap_reduction:UM)
    for (auto f : filenames) {
        // int i = omp_get_thread_num();
        // int n = omp_get_num_threads();
        // printf("hello from thread %d of %d with filename %s\n", i, n, f.c_str());
        compute_file(f, UM);
	}
    // return 0;
    // printMap(UM);
    

    //compute_file(filenames[0], UM);

	auto stop1 = omp_get_wtime();
	
	// sorting in descending order
	ranking rank(UM.begin(), UM.end());

	auto stop2 = omp_get_wtime();
	std::printf("Compute time (s) %f\nSorting time (s) %f\n",
				stop1 - start, stop2 - stop1);
	
	if (showresults) {
		// show the results
		std::cout << "Unique words " << rank.size() << "\n";
		std::cout << "Total words  " << total_words << "\n";
		std::cout << "Top " << topk << " words:\n";
		auto top = rank.begin();
		for (size_t i=0; i < std::clamp(topk, 1ul, rank.size()); ++i)
			std::cout << top->first << '\t' << top++->second << '\n';
	}
}
	
