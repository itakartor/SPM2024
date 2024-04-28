//
// Sequential code of the first SPM Assignment a.a. 23/24.
//
// compile:
// g++ -std=c++20 -O3 -march=native -I<path-to-include> UTWavefront.cpp -o UTW
//
#include <iostream>
#include <vector>
#include <thread>
#include <random>
#include <cassert>
#include <hpc_helpers.hpp>
#include <cmath>

int random(const int &min, const int &max) {
	static std::mt19937 generator(117);
	std::uniform_int_distribution<int> distribution(min,max);
	return distribution(generator);
};		

// emulate some work
void work(std::chrono::microseconds w) {
	auto end = std::chrono::steady_clock::now() + w;
    while(std::chrono::steady_clock::now() < end);	
}

void wavefront(const std::vector<int> &M, const uint64_t &N) {
	for(uint64_t k = 0; k< N; ++k) {        // for each upper diagonal
		
        
        for(uint64_t i = 0; i< (N-k); ++i) {// for each elem. in the diagonal

			work(std::chrono::microseconds(M[i*N+(i+k)])); 
		}
	}
}

void blockCyclicDataDistribution(const uint64_t lineMatrix, const uint64_t& id, 
    const uint64_t num_threads, const uint64_t chunk_size,
    const uint64_t num_tasks,
	std::vector<int>& taskMatrix,
	uint64_t sizeTaskMatrix) {
	
    const uint64_t stride = num_threads*chunk_size;
	const uint64_t offset = id*chunk_size;

	for(uint64_t lower = offset; lower < num_tasks; lower += stride) {
		const uint64_t upper = std::min(lower+chunk_size, num_tasks);
        //std::cout<<"lower: "<<lower<<" upper:"<<upper<<std::endl;
		
		for(uint64_t numElem = lower; numElem < upper; numElem++) {
			//std::cout<<"row: "<<row<<std::endl;
			
			std::cout<<"id: "<<id<<std::endl;
			std::cout<<"numElem: "<<numElem<<std::endl;
			std::cout<<"line: "<<lineMatrix<<std::endl;
			std::cout<<"M: "<<taskMatrix[lineMatrix*sizeTaskMatrix+(lineMatrix+numElem)]<<std::endl;
			//returnVector.emplace_back(row);
		}
		//returnVector.emplace_back(upper-lower);
	}
}

void printMatrix(std::vector<int> &M, int N, int L) {
    for(uint64_t k = 0; k< N; ++k) {  
		for(uint64_t i = 0; i< L; ++i) { 
            std::cout << " [" << k << "," << i << "] " << M[i*N+(i+k)];
    	}
		std::cout << std::endl;
    }
}

int main(int argc, char *argv[]) {
	uint64_t dimMatrix = 5;//512;    // default size of the matrix (NxN)
	int min    = 0;      // default minimum time (in microseconds)
	int max    = 1000;   // default maximum time (in microseconds)
	uint64_t numThreads = 3; // default number of threads
	uint64_t sizeChunck = 2; // default chunks size

	if (argc != 1 && argc != 2 && argc != 4 && argc != 5 && argc != 6) {
		std::printf("use: %s N [min max]\n", argv[0]);
		std::printf("     N size of the square matrix\n");
		std::printf("     min waiting time (us)\n");
		std::printf("     max waiting time (us)\n");		
		std::printf("     num Threads\n");		
		std::printf("     num chunkSize\n");		
		return -1;
	}
	if (argc > 1) {
		dimMatrix = std::stol(argv[1]);
		if (argc > 2) {
			min = std::stol(argv[2]);
			max = std::stol(argv[3]);
		}
		if(argc > 3) {
			numThreads = std::stol(argv[4]);
		}
		if(argc > 4) {
			sizeChunck = std::stol(argv[5]);
		}
	}

	// allocate the matrix
	std::vector<int> M(dimMatrix*dimMatrix, -1);
	// double summation = N*(N+1)/2;
	// std::cout <<"summation: "<< summation <<std::endl;
	
	// double sizeTaskMatrix = round(summation/numThreads);
	
	// std::cout <<"dimTaskMatrix: "<< sizeTaskMatrix <<std::endl;

	// std::vector<int> taskMatrix(numThreads*sizeTaskMatrix, -1);

	uint64_t expected_totaltime=0;
	// init function
	auto init=[&]() {
		for(uint64_t k = 0; k< dimMatrix; ++k) {  
			for(uint64_t i = 0; i< (dimMatrix-k); ++i) {  
				int t = random(min,max);
				M[i*dimMatrix+(i+k)] = t;
				expected_totaltime +=t;				
			}
		}
	};
	
	init();

	std::printf("Estimated compute time ~ %f (ms)\n", expected_totaltime/1000.0);
	
	//TIMERSTART(wavefront);
	
	for(uint64_t numDiagonal = 0; numDiagonal< dimMatrix; ++numDiagonal) {
		for(uint64_t index = 0; index < numThreads; index++) {
			blockCyclicDataDistribution(numDiagonal, index, numThreads, sizeChunck, dimMatrix-numDiagonal, M, dimMatrix);
			std::cout<<"################# <"<<index<<">"<<std::endl;
		}
	}
	printMatrix(M, dimMatrix, dimMatrix);
	//wavefront(M, N); 
    //TIMERSTOP(wavefront);

    return 0;
}
