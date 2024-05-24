// ASSIGNMENT 3

// ----------------------
//
// STRUCTURE:				(files)  (lines)
//	 					|	Source-->Splitter -->Tokenizer--> |	     
//          		    |         			 	  			  |  --> Counter --> | Sink           
//	filelist.txt -->	|	Source-->Splitter -->Tokenizer--> |		                     (map_reduce)	
//                      | 					 	  			  |	 --> Counter --> | Sink
//   					|	Source-->Splitter -->Tokenizer--> |
//  
//     /<------------------------------------ pipeline ------------------------------------->/
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
uint64_t total_words{0};
volatile uint64_t extraworkXline{0};
// ----------------------


//FARM: from filelist.txt to open all files and distribute them among workers
//input:vector of strings		output:#files

//start: filenames.txt --> suddiviso per #nodi (blocchi distribuiti su workers)
//altra idea: dynamic(1,) ma da studiare

struct File_reader: ff_monode_t<std::vector<std::string>, std::vector<std::string>>{ 
	
	const std::vector<std::string> &filenames;
	const uint64_t num_lines;
	volatile uint64_t extraworkXline;
	
	File_reader(const std::vector<std::string> &filenames_, const uint64_t &num_lines_, volatile uint64_t &extraworkXline_) : filenames(filenames_), num_lines(num_lines_), extraworkXline(extraworkXline_) {}

		std::vector<std::string>* svc(std::vector<std::string>*) {

			uint64_t id = get_my_id();
			uint64_t chunk_size = 10;
			uint64_t files = filenames.size();
			
			std::ifstream file(filenames[id], std::ios_base::in); //legge un file alla volta corrispondente all'id del worker
			if (file.is_open()) { //check filename aperto
			std::string line; //lettura
			std::vector<std::string> *V; //crea un vettore di storage per le linee del file
			
			while(!file.eof()){
				V = new std::vector<std::string> ();
				V ->reserve(num_lines);
				uint64_t count = 0;
				while(count<num_lines && std::getline(file, line)) {
					if (!line.empty()) {
						V->emplace_back(line);
						count++;
					}
			}
			// for(auto v = V->begin(); v != V->end(); v++){
			// 	//std::cout << v->c_str() << std::endl; //debug
			// }
			ff_send_out(V); 
			}
		file.close(); 
		return GO_ON;
		}
	return EOS; 
	}
};

//NODE: from lines of a text file to all the tokenized words
//left workers:#lines		right workers:#words per line
struct Word_tokenizer: ff_monode_t<std::vector<std::string>, std::vector<std::string>>{ 

	umap UM;

	Word_tokenizer(umap UM_) : UM(UM_) {}

	std::vector<std::string>* svc(std::vector<std::string> *lines) {

		for(auto line = lines->begin(); line != lines->end(); line++){
			// std::cout << line->c_str() << std::endl; //debug
			char *tmpstr;
			char *token = strtok_r(const_cast<char*>(line->c_str()), " \r\n", &tmpstr);

			while(token) {
				UM[std::string(token)]++;
				token = strtok_r(NULL, " \r\n", &tmpstr);
				total_words++;
			}
			for(volatile uint64_t j {0}; j < extraworkXline; j++);
		}	
		delete lines;
		return GO_ON;
		}

};

//map-reduce: count words' occurences and map them 
// struct Word_counter: ff_monode_t<const std::string&, umap&>{
	
// };


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
				if(filenames.size()<2){
					if (std::filesystem::is_regular_file(line))
					filenames.push_back(line);
				else
					std::cout << line << " is not a regular file, skipt it\n";}
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

	//T
	std::vector<ff_node*> L;

	//temp num_lines (input value?)
	uint64_t num_lines=5;
	// start the time
	ffTime(START_TIME);
	auto start=ffTime(GET_TIME)/1000;

	//call the pipeline that does all that
	// ff_Pipe pipe(File_reader, Word_tokenizer); //, Word_counter);
	// if (pipe.run_and_wait_end()<0) {
    //     error("running pipeMain\n");
    //     return -1;
    // }

	//DAL CODICE DI TORQUATI:

	ff_pipeline* main_pipe = new ff_pipeline(false, qlen, qlen, true);
        
        main_pipe->add_stage(new File_reader(filenames, num_lines, extraworkXline), true);
        main_pipe->add_stage(new Word_tokenizer(UM), true);
        L.push_back(main_pipe);

	if ((*main_pipe).run_and_wait_end()<0) {
        error("running pipeMain\n");
        return -1;
    }
	
	ffTime(STOP_TIME);
	auto stop1 = ffTime(GET_TIME)/1000;
	
	// sorting in descending order
	ranking rank(UM.begin(), UM.end());

	ffTime(STOP_TIME);
	auto stop2 =  ffTime(GET_TIME)/1000;
	
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
	
