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
#include <barrier>
#include <bits/stdc++.h> 

bool showDebug = false;

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

//old wavefront method
void wavefront(const std::vector<int> &M, const uint64_t &N) {
	for(uint64_t k = 0; k< N; ++k) {        // for each upper diagonal
		for(uint64_t i = 0; i< (N-k); ++i) {// for each elem. in the diagonal
			work(std::chrono::microseconds(M[i*N+(i+k)])); 
		}
	}
}

std::vector<uint64_t> blockCyclicDataDistribution(
	const uint64_t lineMatrix, const uint64_t& id, 
    const uint64_t num_threads, const uint64_t chunk_size,
    const uint64_t num_tasks,
	std::vector<int>& taskMatrix,
	uint64_t sizeTaskMatrix) {
	
    const uint64_t stride = num_threads*chunk_size;
	const uint64_t offset = id*chunk_size;

	std::vector<uint64_t> returnedVector;

	for(uint64_t lower = offset; lower < num_tasks; lower += stride) {
		const uint64_t upper = std::min(lower+chunk_size, num_tasks);
        //std::cout<<"lower: "<<lower<<" upper:"<<upper<<std::endl;
		
		for(uint64_t numElem = lower; numElem < upper; numElem++) {
			//std::cout<<"row: "<<row<<std::endl;
			if(showDebug) {
				std::cout<<"id: "<<id<<std::endl;
				std::cout<<"numElem: "<<numElem<<std::endl;
				std::cout<<"line: "<<lineMatrix<<std::endl;
				std::cout<<"M: "<<taskMatrix[lineMatrix*sizeTaskMatrix+(lineMatrix+numElem)]<<std::endl;
			}
			returnedVector.emplace_back(taskMatrix[lineMatrix*sizeTaskMatrix+(lineMatrix+numElem)]);
		}
	}
	returnedVector.emplace_back(-1);

	return returnedVector;
}

void printMatrix(std::vector<int> &M, int N, int L) {
    if(!showDebug) return;
	for(uint64_t k = 0; k< N; ++k) {  
		for(uint64_t i = 0; i< L; ++i) { 
            std::cout << " [" << k << "," << i << "] " << M[i*N+(i+k)];
    	}
		std::cout << std::endl;
    }
}

void printVector(std::vector<uint64_t>& currentVector, uint64_t* test_sum) {
	if(!showDebug) return;
	std::cout << "[ "; 
  
        // Printing vector contents 
        for (auto element : currentVector) {
			if(element != -1) {
				(*test_sum) += element;
			}
            std::cout << element << ' '; 
		}
  
        std::cout << ']'; 
        std::cout << '\n'; 
}  

void printList(std::list<std::vector<uint64_t> >& listOfVectors) { 	
	if(!showDebug) return;
	uint64_t test_sum = 0;
    for (auto vect : listOfVectors) { 
        // Each element of the list is 
        // a vector itself 
        std::vector<uint64_t> currentVector = vect; 
  
        printVector(currentVector, &test_sum);
    }
	std::cout <<"estimate time: " << test_sum << std::endl; 
} 

//use a vector and a common index for read all the task in the task vector
void jobVectorVersion(uint64_t idThread, std::vector<uint64_t> tasks, std::barrier<std::function<void()>>& bar, uint64_t numDiagonal) {
	
	uint64_t test_sum = 0;
	if(showDebug){
		printf("<%ld> print tasks\n", idThread);
		printVector(tasks, &test_sum);
		printf("<%ld> #######################################\n", idThread);
	}
	uint64_t count = 0;
	for(int index = 0; index < tasks.size(); index++) {
		if(tasks[index] == -1) {
			count++;
			if(showDebug)
				printf("<%ld> -1 task, i arrived to the barrier\n", idThread);
			if(count == (numDiagonal-1)) {
				bar.arrive_and_drop();
				return;
			} else {
				bar.arrive_and_wait();
			}
		} else {
			if(showDebug)
				printf("<%ld> i am working on %ld\n", idThread, tasks[index]);
			work(std::chrono::microseconds(tasks[index]));
		}
	}

	if(showDebug)
		printf("<%ld> i finished the tasks, i left\n", idThread);
}

//it's the same version of the vector but use a fifo mode
void jobListVersion(uint64_t idThread, std::vector<uint64_t> tasks, std::barrier<std::function<void()>>& bar, uint64_t numDiagonal) {
	uint64_t test_sum = 0;
	if(showDebug){
		printf("<%ld> print tasks\n", idThread);
		printVector(tasks, &test_sum);
		printf("<%ld> #######################################\n", idThread);
	}
		
	while(!tasks.empty()) {
		if(tasks[0] == -1) {
			tasks.erase(tasks.begin());
			if(showDebug)
				printf("<%ld> -1 task, i arrived to the barrier\n", idThread);
			bar.arrive_and_wait();
		} else {
			if(showDebug)
				printf("<%ld> i am working on %ld\n", idThread, tasks[0]);
			work(std::chrono::microseconds(tasks[0]));
			tasks.erase(tasks.begin());
		}
	}
	printf("<%ld> i finished the tasks, i left\n", idThread);
}

void newWavefront(const uint64_t &dimMatrix, const uint64_t &numThreads, 
	const uint64_t &sizeChunck, std::vector<int>& taskMatrix, std::barrier<std::function<void()>>& myBarrier) {
	std::list<std::vector<uint64_t>> listOfVectorTask;
	std::vector<std::vector<uint64_t>> vectorOfVectorTask;
	std::vector<std::thread> threads;

	for(uint64_t numDiagonal = 0; numDiagonal< dimMatrix; ++numDiagonal) {
		for(uint64_t index = 0; index < numThreads; index++) {
			std::vector<uint64_t> newVector = blockCyclicDataDistribution(numDiagonal, index, numThreads, sizeChunck, dimMatrix-numDiagonal, taskMatrix, dimMatrix);
			if(vectorOfVectorTask.size() < numThreads) {
				vectorOfVectorTask.emplace_back(newVector);
			} else {
				vectorOfVectorTask[index].insert((vectorOfVectorTask[index]).end(), newVector.begin(), newVector.end());
			}

			//this is a list version for merge the task vectors in the main list 
			// if(listOfVectorTask.size() < numThreads) {
			// 	listOfVectorTask.push_back(newVector);				
			// } else {
    		// 	std::list<std::vector<uint64_t>>::iterator it = listOfVectorTask.begin();
			// 	advance(it, index);

			// 	// std::vector<uint64_t> mergeVector((*it).size() + newVector.size());
			// 	// merge((*it).begin(), (*it).end(), newVector.begin(), newVector.end(), mergeVector.begin());
			// 	(*it).insert((*it).end(), newVector.begin(), newVector.end()); 
			// 	//replace(listOfVectorTask.begin(), listOfVectorTask.end(), *it, mergeVector);
			// }
			//std::cout<<"################# <"<<index<<">"<<std::endl;
		}
	}

	if(showDebug){
		printMatrix(taskMatrix, dimMatrix, dimMatrix);
		//std::cout<<"##########################"<<std::endl;
		printList(listOfVectorTask);
	} 

	// for(uint64_t index = 0; index < numThreads; index++) {
	// 	std::list<std::vector<uint64_t>>::iterator it = listOfVectorTask.begin();
	// 	advance(it, index);
	// 	threads.emplace_back(job, index, *it, std::ref(myBarrier));
	// }

	for(uint64_t index = 0; index < numThreads; index++) {
		threads.emplace_back(jobVectorVersion, index, vectorOfVectorTask[index], std::ref(myBarrier), dimMatrix);
	}

	for(auto& thread: threads) {
		thread.join();
	}	
}

int main(int argc, char *argv[]) {
	uint64_t dimMatrix = 512;    // default size of the matrix (NxN)
	int min    = 0;      // default minimum time (in microseconds)
	int max    = 1000;   // default maximum time (in microseconds)
	uint64_t numThreads = std::thread::hardware_concurrency(); // default number of threads
	uint64_t sizeChunck = 2; // default chunks size

	if (argc != 1 && argc != 2 && argc != 4 && argc != 5 && argc != 6 && argc != 7) {
		std::printf("use: %s N [min max] [num Threads] [num chunkSize] [showDebug]\n", argv[0]);
		std::printf("     N size of the square matrix\n");
		std::printf("     min waiting time (us)\n");
		std::printf("     max waiting time (us)\n");		
		std::printf("     num Threads\n");		
		std::printf("     num chunkSize\n");		
		std::printf("     prints of showDebug\n");		
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
		if(argc > 5) {
			int tmp;
			try { 
				tmp=std::stol(argv[6]);
			} catch(std::invalid_argument const& ex) {
				std::printf("%s is an invalid number (%s)\n", argv[6], ex.what());
				return -1;
			}
			if (tmp == 1) showDebug = true;
		}
	}

	std::barrier<std::function<void()>> myBarrier(numThreads, []() noexcept 
     {
		if(showDebug)
			std::cout << "barrier reached"<<std::endl;
	});

	// allocate the matrix
	std::vector<int> M(dimMatrix*dimMatrix, -1);
	
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
	
	TIMERSTART(wavefront);
	wavefront(M, dimMatrix); 
    TIMERSTOP(wavefront);

	TIMERSTART(newWavefront);
	newWavefront(dimMatrix, numThreads, sizeChunck, M, std::ref(myBarrier)); 
    TIMERSTOP(newWavefront);

}
