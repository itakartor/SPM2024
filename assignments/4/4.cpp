#include <cstdio>
#include <random>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include "mpi.h"

//COMPILE
//$ mpicxx -std=c++20 parallel_nkeys.cpp -o nkeys

const long   SIZE = 64;

struct params{
	int m, k, n;
	float alpha;
};

void readParams(params* p, const long c1, const long c2){
	// reading here from a file
	p->m = c1;
	p->k = c2;
	p->n = c1;
	p->alpha = 1;
}



long random(const int &min, const int &max) {
	static std::mt19937 generator(117);
	std::uniform_int_distribution<long> distribution(min,max);
	return distribution(generator);
};		


void init(float *M, const long c1, const long c2, const long key) {
	// for(long i=0;i<c1;++i)
	// 	for(long j=0;j<c2;++j)
	// 		M[i][j] = (key-i-j)/static_cast<float>(SIZE);
	for(long i=0; i<c1; i++)
        for(long j=0; j<c2; j++)
            M[i*c2+j] = (key-i-j)/static_cast<float>(SIZE);
}

// matrix multiplication:  C = A x B  A[c1][c2] B[c2][c1] C[c1][c1]
// mm returns the sum of the elements of the C matrix
auto mm(const auto& myA, const auto& myB, int blockRows, int blockCols, auto& params) {
	float sum{0};
	for(int i=0; i<blockRows; i++){
		for(int j=0; j<blockCols; j++){
			// sum[i*blockCols+j] = 0.0;
			auto accum = float(0.0);
			for(int l=0; l<params.k; l++){
				// myC[i*blockCols+j]
				accum += params.alpha*myA[i*params.k+l]*myB[l*blockCols+j];
			}	
		sum += accum;
		}
	}
	return sum;
}

// initialize two matrices with the computed values of the keys
// and execute a matrix multiplication between the two matrices
// to obtain the sum of the elements of the result matrix 
float compute(const long c1, const long c2, int gridDim, long key1, long key2, int myId, MPI_Datatype paramsType) {
	params p;

	if(!myId){
		// Process 0 reads the parameters from a configuration file
		readParams(&p, c1, c2);
	}
	// Broadcast of all the parameters using one message with a struct
	MPI_Bcast(&p, 1, paramsType, 0, MPI_COMM_WORLD);

	if((p.m < 1) || (p.n < 1) || (p.k<1)){
		// Only the first process prints the output message
		if(!myId){
			std::cout << "ERROR: 'm', 'k' and 'n' must be higher than 0" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, -1);
		}	
	}

	if((p.m%gridDim) || (p.n%gridDim)){
		// Only the first process prints the output message
		if(!myId){
			std::cout << "ERROR: 'm', 'n' must be multiple of the grid dimensions" << std::endl;
		}
		MPI_Abort(MPI_COMM_WORLD, -1);
	}


	float *A=nullptr;
	float *B=nullptr;
	// float *C=nullptr;
	float *myA;
	float *myB;
	// float *myC;

	float tot_sum{0};
	int blockRows = p.m/gridDim;
	int blockCols = p.n/gridDim;
	MPI_Request req;
	
	// Only one process reads the data from the files
	if(!myId){
		A = new float[p.m*p.k];
		B = new float[p.k*p.n];
		// std::vector<std::vector<float>> A(c1, std::vector<float>(c2,0.0));
		// std::vector<std::vector<float>> B(c2, std::vector<float>(c1,0.0));
		init(A, c1, c2, key1);
		init(B, c2, c1, key2);
		// readInput(p.m, p.k, p.n, A, B);
	}
	MPI_Barrier(MPI_COMM_WORLD); //every process has A and B



	// Create the datatype for a block of rows of A
	MPI_Datatype rowsType;
	MPI_Type_contiguous(blockRows*p.k, MPI_FLOAT, &rowsType);
	MPI_Type_commit(&rowsType);

	// Send the rows of A that needs each process
	if(!myId){
		for(int i=0; i<gridDim; i++)
			for(int j=0; j<gridDim; j++)
				MPI_Isend(&A[i*blockRows*p.k], 1, rowsType, i*gridDim+j, 0, MPI_COMM_WORLD, &req);
	}	
	myA = new float[blockRows*p.k];
	MPI_Recv(myA, 1, rowsType, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);				

	// Create the datatype for a block of columns of B
	MPI_Datatype colsType;
	MPI_Type_vector(p.k, blockCols, p.n, MPI_FLOAT, &colsType);
	MPI_Type_commit(&colsType);

	// Send the columns of B that needs each process
	if(!myId){
		for(int i=0; i<gridDim; i++)
			for(int j=0; j<gridDim; j++)
				MPI_Isend(&B[blockCols*j], 1, colsType, i*gridDim+j, 0, MPI_COMM_WORLD, &req);
	}
		
	myB = new float[p.k*blockCols];
	MPI_Recv(myB, p.k*blockCols, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	tot_sum += mm(myA, myB, blockRows, blockCols, p);


	if(!myId){
		float tmp{0};
		// for(int i=0; i<blockRows; i++)
		// 	memcpy(&C[i*p.n], &myC[i*blockCols], blockCols*sizeof(float));		
		for(int i=0; i<gridDim; i++)
			for(int j=0; j<gridDim; j++)
				if(i || j){
					MPI_Recv(&tmp, 1, MPI_FLOAT, i*gridDim+j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					tot_sum += tmp;
				}
	} else 
		MPI_Send(&tot_sum, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

	//free types
	// Delete the types and arrays
	MPI_Type_free(&rowsType);
	MPI_Type_free(&colsType);

	if(!myId){
		return tot_sum;
	}
	else{
		return 0;
	}
} 




int main(int argc, char* argv[]) {
	if (argc < 3) {
		std::printf("use: %s nkeys length [print(0|1)]\n", argv[0]);
		std::printf("     print: 0 disabled, 1 enabled\n");
		return -1;
	}
	//INITIALIZATION
	MPI_Init(&argc, &argv);

	int myId;
	int numP;
	// Get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &myId); 
    MPI_Comm_size(MPI_COMM_WORLD, &numP); // == length? 
	int gridDim = sqrt(numP);

	if(gridDim*gridDim != numP){
		// Only the first process prints the output message
		if(!myId)
			std::cout << "ERROR: the number of processes must be square" << std::endl;

		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	

	long nkeys  = std::stol(argv[1]);  // total number of keys
	// length is the "stream length", i.e., the number of random key pairs
	// generated
	long length = std::stol(argv[2]);  
	bool print = false;
	if (argc == 4)
		print = (std::stoi(argv[3]) == 1) ? true : false;
	
	long key1, key2;

	std::map<long, long> map;
	for(long i=0;i<nkeys; ++i) map[i]=0;
	
	std::vector<float> V(nkeys, 0);
	bool resetkey1=false;
	bool resetkey2=false;


	// DATATYPES
	int blockLengths[2] = {3, 1};
	MPI_Aint lb, extent;
	MPI_Type_get_extent(MPI_INT, &lb, &extent); 
	MPI_Aint disp[2] = {0, 3*extent};
	MPI_Datatype types[2] = {MPI_INT, MPI_FLOAT}; 

	// Create the datatype for the parameters
	MPI_Datatype paramsType;
	MPI_Type_create_struct(2, blockLengths, disp, types, &paramsType);
	MPI_Type_commit(&paramsType);


	//START MAP FILLING
	for(int i=0;i<length; ++i) {
		key1 = random(0, nkeys-1);  // value in [0,nkeys[
		key2 = random(0, nkeys-1);  // value in [0,nkeys[
		
		
		if (key1 == key2) // only distinct values in the pair
			key1 = (key1+1) % nkeys; 

		std::cout << "TEST PRINT: \n key 1: " << key1 << "key 2: " << key2 << std::endl;

		map[key1]++;  // count the number of key1 keys
		map[key2]++;  // count the number of key2 keys

		float r1;
		float r2;
		// if key1 reaches the SIZE limit, then do the computation and then
		// reset the counter ....
		if (map[key1] == SIZE && map[key2]!=0) {			
			r1= compute(map[key1], map[key2], gridDim, key1, key2, myId, paramsType);
			V[key1] += r1;  // sum the partial values for key1
			resetkey1=true;			
		}
		// if key2 reaches the SIZE limit ....
		if (map[key2] == SIZE && map[key1]!=0) {			
			r2= compute(map[key2], map[key1], gridDim, key2, key1, myId, paramsType);
			V[key2] += r2;  // sum the partial values for key1
			resetkey2=true;
		}
		if (resetkey1) {
			// updating the map[key1] initial value before restarting
			// the computation
			auto _r1 = static_cast<unsigned long>(r1) % SIZE;
			map[key1] = (_r1>(SIZE/2)) ? 0 : _r1;
			resetkey1 = false;
		}
		if (resetkey2) {
			// updating the map[key2] initial value before restarting
			// the computation
			auto _r2 = static_cast<unsigned long>(r2) % SIZE;
			map[key2] = (_r2>(SIZE/2)) ? 0 : _r2;
			resetkey2 = false;
		}
	}

	// compute the last values
	for(long i=0;i<nkeys; ++i) {
		for(long j=0;j<nkeys; ++j) {
			if (i==j) continue;
			if (map[i]>0 && map[j]>0) {
				auto r1= compute(map[i], map[j], gridDim, i, j, myId, paramsType);
				auto r2= compute(map[j], map[i], gridDim, j, i, myId, paramsType);
				V[i] += r1;
				V[j] += r2;
				std::cout << "TEST PRINT: \n r1: " << r1 << " r2: " << r2 << std::endl;

			}
		}
	}
	// printing the results
	if (print) {
		for(long i=0;i<nkeys; ++i)
			std::printf("key %ld : %f\n", i, V[i]);
	}
	
	MPI_Type_free(&paramsType);
	// Terminate MPI
	MPI_Finalize();
}