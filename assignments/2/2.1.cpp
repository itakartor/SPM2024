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

// compile: g++ -std=c++17 -O3 -march=native -fopenmp 2.1.cpp -o 2.1

using umap=std::unordered_map<std::string, uint64_t>;
using pair=std::pair<std::string, uint64_t>;
struct Comp {
	bool operator ()(const pair& p1, const pair& p2) const {
		return p1.second > p2.second;
	}
};
using ranking=std::multiset<pair, Comp>;

// ------ globals --------
uint64_t total_words{0};
volatile uint64_t extraworkXline{0};
uint64_t maxNumLines = 20;
// ----------------------

void printMap(umap myMap) {
    printf("MAP:\n");
    for(auto& elem : myMap) {
        printf("%s, %ld\n", elem.first.c_str(), elem.second);
    }
}

void tokenize_line2(const std::string& line, umap& UM) {
	char *tmpstr;
	char *token = strtok_r(const_cast<char*>(line.c_str()), " \r\n", &tmpstr);
	while(token) {
		++UM[std::string(token)];
		token = strtok_r(NULL, " \r\n", &tmpstr);
		
		#pragma omp atomic
		++total_words;
	}
	//for(volatile uint64_t j{0}; j<extraworkXline; j++);
}

void tokenize_line(const std::string& line, std::vector<umap>& UM) {
	char *tmpstr;
	char *token = strtok_r(const_cast<char*>(line.c_str()), " \r\n", &tmpstr);
	while(token) {
		++UM[omp_get_thread_num()][std::string(token)];
		token = strtok_r(NULL, " \r\n", &tmpstr);
		
		#pragma omp atomic
		++total_words;
	}
	//for(volatile uint64_t j{0}; j<extraworkXline; j++);
}

void compute_file(const std::string& filename, std::vector<umap>& UMVector) {
	std::ifstream file(filename, std::ios_base::in);

	if (file.is_open()) {
		std::string readingLine;

		while(!file.eof()) {
			std::vector<std::string>* V = new std::vector<std::string>();
			V->reserve(maxNumLines);
			uint64_t counter = 0;

			while(counter < maxNumLines && !(std::getline(file, readingLine).eof())) {
				if (!readingLine.empty()) {
					V->emplace_back(readingLine);
					counter++;
				}
			}

			#pragma omp task firstprivate(V) shared(UMVector)
			{	
				// int i = omp_get_thread_num();
				// int n = omp_get_num_threads();
				// printf("task %d/%d\n", i,n);
				// umap UM;
				for(auto line=V->begin(); line != V->end(); line++) {
					// tokenize_line(*line, UM);
					tokenize_line(*line, UMVector);
				}
				// printMap(UM);
					
				delete V;

				// #pragma omp critical
				// UMVector.emplace_back(UM);
			}
		}
		

	} 
	file.close();
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

int main(int argc, char *argv[]) {

	auto usage_and_exit = [argv]() {
		std::printf("use: %s filelist.txt [extraworkXline] [topk] [showresults] [maxNumLines]\n", argv[0]);
		std::printf("     filelist.txt contains one txt filename per line\n");
		std::printf("     extraworkXline is the extra work done for each line, it is an integer value whose default is 0\n");
		std::printf("     topk is an integer number, its default value is 10 (top 10 words)\n");
		std::printf("     showresults is 0 or 1, if 1 the output is shown on the standard output\n");
		std::printf("     maxNumLines is the max number of line in a single block for the tokenize task. the default value is 20\n\n");
		exit(-1);
	};

	std::vector<std::string> filenames;
	size_t topk = 10;
	bool showresults = false;
	if (argc < 2 || argc > 6) {
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
			if (argc >= 5) {
				int tmp;
				try { 
					tmp=std::stol(argv[4]);
				} catch(std::invalid_argument const& ex) {
					std::printf("%s is an invalid number (%s)\n", argv[4], ex.what());
					return -1;
				}
				if (tmp == 1) showresults = true;
			}
			if(argc == 6) {
				try { 
					maxNumLines = std::stoul(argv[5]);
				} catch(std::invalid_argument const& ex) {
					std::printf("%s is an invalid number (%s)\n", argv[5], ex.what());
					return -1;
				}
				if (maxNumLines == 0) {
					std::printf("%s must be a positive integer\n", argv[5]);
					return -1;
				}
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

	// start the time
	auto start = omp_get_wtime();
	// uint64_t num_threads = omp_get_max_threads();	
	// std::cout << "max num threads: " << num_threads << std::endl;

	// used for storing results
	std::vector<umap> vectorUM(extraworkXline);
	umap UM;
	// std::cout<<"num lines: " << maxNumLines << std::endl;
	#pragma omp parallel for schedule(dynamic, 1) num_threads(extraworkXline) 
	for (auto f : filenames) {
		// int i = omp_get_thread_num();
		// int n = omp_get_num_threads();
		// printf("file %d/%d\n", i,n);
		compute_file(f, vectorUM);
	}
	
	for(auto& i : vectorUM) {
		//printMap(i);
		reduce_umaps(UM, i);
	}

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
	
