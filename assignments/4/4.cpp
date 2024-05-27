#include <cstdio>
#include <random>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include "mpi.h"

//COMPILE
//$ mpicxx -Wall -O3 -std=c++20 parallel_nkeys.cpp -o nkeys
//ALLOCATE
//$ salloc -N 4 --time=00:30:00
//SLURM
//$ mpirun -n 4 --host node01,node02,node03,node04 ./nkeys 2 100 1
//$ srun --mpi=pmix -n 4 ./newnkeys 2 12 1

const long   SIZE = 64;

bool debug=false;

struct params{
	int m, k, n;
	float alpha;
};

struct keypair{
	float key1, key2;
};

struct computation{
	float m_1, m_2, key1, key2;
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


void init(auto& M, const long c1, const long c2, const long key) {
	for(long i=0;i<c1;++i)
		for(long j=0;j<c2;++j)
			M[i][j] = (key-i-j)/static_cast<float>(SIZE);
}

// matrix multiplication:  C = A x B  A[c1][c2] B[c2][c1] C[c1][c1]
// mm returns the sum of the elements of the C matrix
auto mm(const auto& A, const auto& B, const long c1,const long c2) {
	float sum{0};
	
    for (long i = 0; i < c1; i++) {
        for (long j = 0; j < c1; j++) {
            auto accum = float(0.0);
            for (long k = 0; k < c2; k++)
                accum += A[i][k] * B[k][j];
            sum += accum;
		}
	}
	return sum;
}

// initialize two matrices with the computed values of the keys
// and execute a matrix multiplication between the two matrices
// to obtain the sum of the elements of the result matrix 
float compute(const long c1, const long c2, long key1, long key2) {
	std::vector<std::vector<float>> A(c1, std::vector<float>(c2,0.0));
	std::vector<std::vector<float>> B(c2, std::vector<float>(c1,0.0));

	init(A, c1, c2, key1);
	init(B, c2, c1, key2);
	auto r = mm(A,B, c1,c2);

	return r;
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


	long nkeys  = std::stol(argv[1]);  // total number of keys
	// length is the "stream length", i.e., the number of random key pairs generated
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

	//KEYS DATATYPE
	int blockLengths_k[1] = {2};
	MPI_Aint lb_k, extent_k;
	MPI_Type_get_extent(MPI_FLOAT, &lb_k, &extent_k); 
	MPI_Aint disp_k[2] = {0, 2*extent_k};
	MPI_Datatype types_k[1] = {MPI_FLOAT}; 

	// Create the datatype for the parameters
	MPI_Datatype keysType;
	MPI_Type_create_struct(1, blockLengths_k, disp_k, types_k, &keysType);
	MPI_Type_commit(&keysType);

	//COMP DATATYPE
	int blockLengths_c[1] = {4};
	MPI_Aint lb_c, extent_c;
	MPI_Type_get_extent(MPI_FLOAT, &lb_c, &extent_c); 
	MPI_Aint disp_c[2] = {0, 2*extent_c};
	MPI_Datatype types_c[1] = {MPI_FLOAT}; 

	// Create the datatype for the parameters
	MPI_Datatype compType;
	MPI_Type_create_struct(1, blockLengths_c, disp_c, types_c, &compType);
	MPI_Type_commit(&compType);

	MPI_Datatype rowType;
	MPI_Type_contiguous(nkeys, MPI_FLOAT, &rowType);
	MPI_Type_commit(&rowType);

	keypair k;
	computation c;

	MPI_Request req;


	double start = MPI_Wtime();

	//START MAP FILLING
	
	if((length%(numP-1))!=0){
		if(!myId)
			std::cout << "ERROR: length must be a multiple of the number of processes" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	int chunkSize=length/(numP-1);
	if(myId){
		for(int i=0;i<chunkSize; ++i) {
			key1 = random(0, nkeys-1);  // value in [0,nkeys[
			key2 = random(0, nkeys-1);  // value in [0,nkeys[
			if (key1 == key2) // only distinct values in the pair
				key1 = (key1+1) % nkeys; 

			k.key1=key1;
			k.key2=key2;

			MPI_Isend(&k, 1, keysType, 0, 100, MPI_COMM_WORLD, &req);
			if(debug){
				std::cout << "process number: " << myId << "\t keys: " << k.key1 << k.key2 << "\t i= " << i << std::endl;
			}
		} 
		if(debug){
			std::cout << "END key generation for node " << myId << std::endl;
		}
	} 
	else {
		for(int i=0;i<length; ++i){	
			MPI_Recv(&k, 1, keysType, MPI_ANY_SOURCE, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			key1=k.key1;
			key2=k.key2;
			if(debug){
				std::cout << "Received keys: " << key1 << key2 << std::endl;
			}
			map[key1]++;  // count the number of key1 keys
			map[key2]++;  // count the number of key2 keys


			float r1;
			float r2;
			// if key1 reaches the SIZE limit, then do the computation and then
			// reset the counter ....
			if (map[key1] == SIZE && map[key2]!=0) {			
				r1= compute(map[key1], map[key2], key1, key2);
				V[key1] += r1;  // sum the partial values for key1
				resetkey1=true;
				if(debug){
					std::cout << "process " << myId << "\t r1: " << r1 << std::endl;	
				}
			}
			// if key2 reaches the SIZE limit ....
			if (map[key2] == SIZE && map[key1]!=0) {			
				r2= compute(map[key2], map[key1], key2, key1);
				V[key2] += r2;  // sum the partial values for key1
				resetkey2=true;
				if(debug){
					std::cout << "process " << myId << "\t r2: " << r2 << std::endl;
				}
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
		if(debug){
			std::cout << "END key sharing and compute" << std::endl;
		}
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// for(long i=0; i<nkeys; ++i){
	// 		std::cout << "map: " << map[i] << "\t i: "<< i << std::endl;
	// }

	// compute the last values
	if(!myId){
		int roundrobin=1; //id a cui stiamo mandando il messaggio
		if(debug){
			std::cout << "sono nel server" << std::endl;
		}
		for(long i=0;i<nkeys; ++i) {
			for(long j=0;j<nkeys; ++j) {
				if (i==j) continue;
				if (map[i]>0 && map[j]>0) {
					if(roundrobin>(numP-1)){
						roundrobin=1; 
						if(debug){
							std::cout << "roundrobin resettato a 1" << std::endl;
						}
					}

					c.m_1=map[i];
					c.m_2=map[j];
					c.key1=i;
					c.key2=j;

				if(debug){
					std::cout << "PRINT OF c1: key1 = " << c.key1 << " key2= " << c.key2 << " m1= " << c.m_1 << " m2= " << c.m_2 << std::endl;
				}

					MPI_Isend(&c, 1, compType, roundrobin, 101, MPI_COMM_WORLD, &req);
					++roundrobin;

					if(roundrobin>(numP-1)){
						roundrobin=1; 
					}
					
					c.m_1=map[j];
					c.m_2=map[i];
					c.key1=j;
					c.key2=i;

				if(debug){
					std::cout << "PRINT OF c2: key1 = " << c.key1 << " key2= " << c.key2 << " m1= " << c.m_1 << " m2= " << c.m_2 << std::endl;
				}

					MPI_Isend(&c, 1, compType, roundrobin, 101, MPI_COMM_WORLD, &req);
					++roundrobin;
				}
			}
		} 

		c.key1=-1;
		for(int g=1; g<(numP); g++){
			if(debug){
				std::cout << "sono nodo " << myId << " e sto inviando " <<c.key1 << " a nodo " << g << std::endl;
			}
			MPI_Isend(&c, 1, compType, g, 101, MPI_COMM_WORLD, &req);	
		}
		
		std::vector<float> tempV(nkeys, 0);
		for(int j=0; j<(numP-1); j++){
			MPI_Recv(tempV.data(), nkeys, MPI_FLOAT, MPI_ANY_SOURCE, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for(int i=0; i<nkeys;i++){
				if(debug){
					std::cout << "V[i]: " << V[i] <<" temp: "<<tempV[i] << std::endl;
				}

				V[i]+=tempV[i];
				if(debug){
					std::cout << "V[i]: " << V[i] << std::endl;
				}
			}
		}
	} else { //we leave the server and enter in the nodes
		if(debug){
			std::cout << "sono nel nodo " << myId << std::endl;
		}
		bool check=true;
		while(check){
			if(debug){
				std::cout << "sono nodo " << myId << " e sto aspettando di ricevere" << std::endl;
			}
			MPI_Recv(&c, 1, compType, 0, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if(debug){
				printf("sono nodo %d e ho ricevuto: %f, %f, %f, %f", myId, c.key1, c.key2, c.m_1, c.m_2);
			}
			if(c.key1==-1){
				check=false;
				if(debug){
					std::cout << myId << " sta uscendo" << std::endl;
				}
			} else{
				auto r = compute(c.m_1, c.m_2, c.key1, c.key2);
				V[c.key1] += r;
				if(debug){
					std::cout << "V[key1] " << V[key1] << " with key: " << key1 << "\t r: "<< r << std::endl;	
				}
			}	
		}
		if(debug){
			std::cout << "sono nodo " << myId << " e sto mandando V" << std::endl;
		}
		MPI_Isend(V.data(), V.size(), MPI_FLOAT, 0, 102, MPI_COMM_WORLD, &req);
		
		if(debug){
			for(int h=0; h<nkeys; h++){
				std::cout << "V[h]: " << V[h] << std::endl;
			}
		}
	}
	

	double end = MPI_Wtime();

	// printing the results
	if (print && !myId) {
		for(long i=0;i<nkeys; ++i)
			std::printf("Node printing : %d \n key %ld : %f\n", myId, i, V[i]);
	}
	std::cout << "Time with " << numP << " processes: " << end-start << " seconds" << std::endl;
	
	MPI_Type_free(&paramsType);
	MPI_Type_free(&keysType);
	MPI_Type_free(&compType);
	MPI_Type_free(&rowType);

	// Terminate MPI
	MPI_Finalize();
	return 0;

}

//MPI TAGS
// 100: send key generator
// 101: compute data (with also roundrobin)
// 102: vector V data from each node to be reduced/merged


// values: 2 120 1 WITH 4 NODES (1st test)
//PARALLEL TIME: ~2,9 ms
//SEQUENTIAL TIME: ~14 ms

//values: 2 1.000.000 
// SEQUENTIAL TIME: ~144476 ms
// WITH 11 NODES
// PARALLEL TIME: ~27.4509 s
// WITH 21 NODES
// PARALLEL TIME: ~16.9323 s

//values: 100 1.000.000 
// SEQUENTIAL TIME: ~86749 ms
// WITH 11 NODES
// PARALLEL TIME: ~ 27.4736 s
// WITH 21 NODES
// PARALLEL TIME: ~11.9632 s
//WITH 51 NODES 
//PARALLEL TIME: ~7.9109 s