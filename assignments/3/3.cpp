// ASSIGNMENT 3

//$ g++ -I ./fastflow -DNO_DEFAULT_MAPPING  3.cpp -o 3

// ----------------------
//
// STRUCTURE:				(files)  (lines) AllToAll
//	 					|	Reader --> Tokenizer --> |   
//          		    |         			 	  	 |		             
//	filelist.txt -->	|	Reader --> Tokenizer --> |	 reducer in master node                    	
//                      | 					 	  	 |	(map_reduce)	  
//   					|	Reader --> Tokenizer --> |
//  
//	
// ----------------------

#if !defined(DEFAULT_BUFFER_CAPACITY)
#define DEFAULT_BUFFER_CAPACITY 2048
#endif

#include <omp.h> 
#include <cstring>
#include <vector>
#include <set>
#include <string>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <ff/ff.hpp>
#include <ff/pipeline.hpp>

const size_t qlen = DEFAULT_BUFFER_CAPACITY;

using namespace ff;

using umap=std::unordered_map<std::string, uint64_t>;
using pair=std::pair<std::string, uint64_t>;
struct Comp {
	bool operator ()(const pair& p1, const pair& p2) const {
		return p1.second > p2.second;
	}
};
using ranking=std::multiset<pair, Comp>;

// ------ globals --------
volatile uint64_t extraworkXline{0};
volatile uint64_t num_lines{0};
// ----------------------

//struct of the readers that takes the its file list as parameter and the max num lines to sent to the tokenizers nodes 
struct File_reader: ff_monode_t<std::vector<std::string>>{ 
	
	const std::vector<std::string> filenames;
	const uint64_t num_lines;
	File_reader(const std::vector<std::string> filenames_, const uint64_t num_lines_) : filenames(filenames_), num_lines(num_lines_) {}

		std::vector<std::string>* svc(std::vector<std::string>*) {

			// uint64_t id = get_my_id();
			uint64_t numFiles = filenames.size();
			// std::cout << "id: " << id << " filenames " << files << std::endl;
			
			for(uint64_t index = 0 ; index<numFiles;index++){
				// std::cout << "indice filename: "<<filenames[index] << std::endl;
				std::ifstream file(filenames[index], std::ios_base::in);
				if (file.is_open()) { //check filename aperto
					std::string line;
					std::vector<std::string> *V;
					
					while(!file.eof()){
						V = new std::vector<std::string> ();
						V ->reserve(num_lines);
						uint64_t count = 0;
						while(count<num_lines && !std::getline(file, line).eof()) {
							if (!line.empty()) {
								V->emplace_back(line);
								count++;
							}
						}

						ff_send_out(V); 
					}
					file.close(); 
			//return GO_ON;
				}
			}
	return EOS;
	// return GO_ON;  
	}
};

// NODE: from lines of a text file to all the tokenized words
// left workers:#lines		right workers:counts the words
// struct of the tokenizers that takes the its line as parameter and counts and save occurency the worlds in the local umap 
struct Word_tokenizer: ff_node_t<std::vector<std::string>, umap>{ 

	umap UM;
	uint64_t loc_total_words{0};

	Word_tokenizer() {}

	umap*	svc(std::vector<std::string> *lines) {

		// std::cout << (*lines).size() <<std::endl;
		for(auto line = lines->begin(); line != lines->end(); line++){
		// 	// std::cout << line->c_str() << std::endl; //debug
			
			char *tmpstr;
			char *token = strtok_r(const_cast<char*>(line->c_str()), " \r\n", &tmpstr);

			while(token) {
				// std::cout << token << std::endl;
				UM[std::string(token)]++;
				token = strtok_r(NULL, " \r\n", &tmpstr);
				loc_total_words++;
			}
		}	
		// ff_send_out(UM);
		delete lines;
		return GO_ON;
	}
};

// it's a function for reduce two map in only one called output
void reduce_umaps(umap& output, umap& input)
{
	if(!output.empty() || !input.empty()){
		for(auto& X : input) {
			auto t = output.find(X.first);
			if (t == output.end()) {
				output.insert(pair{X.first, X.second});
			} else {
				output.at(X.first) += X.second; //Will throw if X.first doesn't exist in output. 
			}
		}
	}
}

// method for split the tasks from the threads
// it's static method for spitting the tasks if you know the total number of them
// returns a vector of tasks for a single thread and puts -1 
// for limit the end of the tasks of the specific diagonal
std::vector<std::string> blockCyclicDataDistribution(
	const uint64_t id, 
    const uint64_t num_threads, //#worker tot
	const uint64_t chunk_size,
    const uint64_t num_files,
	std::vector<std::string>& fileVector) {
	
    const uint64_t stride = num_threads*chunk_size;
	const uint64_t offset = id*chunk_size;

	std::vector<std::string> returnedVector;

	for(uint64_t lower = offset; lower < num_files; lower += stride) {
		const uint64_t upper = std::min(lower+chunk_size, num_files);
        //std::cout<<"lower: "<<lower<<" upper:"<<upper<<std::endl;
		
		for(uint64_t numElem = lower; numElem < upper; numElem++) {
			returnedVector.emplace_back(fileVector[numElem]);
		}
	}

	return returnedVector;
}


//DISABILITA MAPPING NATIVO FF
//taskset -c   OPENMP
//mapping_string -d FF
int main(int argc, char *argv[]) {

	auto usage_and_exit = [argv]() {
		std::printf("use: %s filelist.txt [extraworkXline] [numlines] [topk] [showresults]\n", argv[0]);
		std::printf("     filelist.txt contains one txt filename per line\n");
		std::printf("     extraworkXline is the extra work done for each line, it is an integer value whose default is 0\n");
		std::printf("     topk is an integer number, its default value is 10 (top 10 words)\n");
		std::printf("     showresults is 0 or 1, if 1 the output is shown on the standard output\n\n");
		exit(-1);
	};
	
	uint64_t num_cores = ff_numCores();
	std::vector<std::string> filenames;
	size_t topk = 10;
	bool showresults=false;


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
			try {num_lines=std::stoul(argv[3]);
			} catch(std::invalid_argument const& ex) {
				std::printf("%s is an invalid number (%s)\n", argv[3], ex.what());
				return -1;
		}
		if (argc > 4) {
			try { topk=std::stoul(argv[4]);
			} catch(std::invalid_argument const& ex) {
				std::printf("%s is an invalid number (%s)\n", argv[4], ex.what());
				return -1;
			}
			if (topk==0) {
				std::printf("%s must be a positive integer\n", argv[4]);
				return -1;
			}
			if (argc == 6) {
				int tmp;
				try { tmp=std::stol(argv[5]);
				} catch(std::invalid_argument const& ex) {
					std::printf("%s is an invalid number (%s)\n", argv[5], ex.what());
					return -1;
				}
				if (tmp == 1) showresults = true;
			}
		}
	}
	}
	
	if (std::filesystem::is_regular_file(argv[1])) {
		std::ifstream file(argv[1], std::ios_base::in);
		if (file.is_open()) {
			std::string line;
			while(std::getline(file, line)) {
				// if(filenames.size()<2){
				if (std::filesystem::is_regular_file(line))
					filenames.push_back(line);
				else
					std::cout << line << " is not a regular file, skip it\n";}
				// }					
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
	ffTime(START_TIME);

	ff_a2a a2a;
	std::vector<ff_node*> left_w;
    std::vector<ff_node*> right_w;

	// std::cout << filenames.size() << std::endl;
	uint64_t chunks=10;
	for (uint64_t i = 0; i < extraworkXline; i++){
		// std::cout << "File_dist: " << files_dist.size() <<std::endl;
		left_w.push_back(new File_reader(blockCyclicDataDistribution(i, extraworkXline, chunks, filenames.size(), filenames), num_lines));
	}
		

	for (uint64_t i = 0; i < num_cores; i++) {
		right_w.push_back(new Word_tokenizer());
	}

    a2a.add_firstset(left_w, 1);
    a2a.add_secondset(right_w);
    
    if (a2a.run_and_wait_end() < 0) {
		error("running a2a\n");
		return -1;
    }

	// defined umap and total words variable
	umap UM_def = ((Word_tokenizer*) right_w[0])->UM;
	uint64_t total_words_main = ((Word_tokenizer*) right_w[0])->loc_total_words;

	// reduce the results in only two main variable
	for(uint64_t i=1; i<right_w.size();i++){
		reduce_umaps(UM_def, (((Word_tokenizer*) right_w[i])->UM));
		total_words_main += ((Word_tokenizer*) right_w[i])->loc_total_words;
	}
	

	
	ffTime(STOP_TIME);
	auto stop1 = ffTime(GET_TIME)/1000;
	
	// sorting in descending order
	ranking rank(UM_def.begin(), UM_def.end());

	ffTime(STOP_TIME);
	auto stop2 =  ffTime(GET_TIME)/1000;
	
	std::printf("Compute time (s) %f\nSorting time (s) %f\n", stop1, stop2 - stop1);
	
	if (showresults) {
		// show the results
		std::cout << "Unique words " << rank.size() << "\n";
		std::cout << "Total words  " << total_words_main << "\n";
		std::cout << "Top " << topk << " words:\n";
		auto top = rank.begin();
		for (size_t i=0; i < std::clamp(topk, 1ul, rank.size()); ++i)
			std::cout << top->first << '\t' << top++->second << '\n';
	}
}
	