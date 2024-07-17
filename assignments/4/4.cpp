#include <cstdio>
#include <random>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "mpi.h"

//MPI TAGS
// 100: send key generator
// 101: compute data (with also roundrobin)
// 102: vector V data from each node to be reduced/merged
#define KEY_GENERATOR 100
#define COMPUTE_DATA 101
#define SEND_VECTOR_V 102

//COMPILE
//$ mpicxx -Wall -O3 -std=c++20 parallel_nkeys.cpp -o nkeys
//ALLOCATE
//$ salloc -N 4 --time=00:30:00
//SLURM
//$ mpirun -n 4 --host node01,node02,node03,node04 ./nkeys 2 100 1
//$ srun --mpi=pmix -n 4 ./newnkeys 2 12 1

const long   SIZE = 64;

bool debug = false;
bool debug2 = false;

//the struct used for send from generation nodes to server 
struct keypair{
	long key1, key2;
};

//the struct used for send information to nodes for run the compute method
struct computation{
	long m_1, m_2, key1, key2;
	float result;
};

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

	std::vector<std::vector<float>> A(c1, std::vector<float>(c2,0.0)); // c1 * c2
	std::vector<std::vector<float>> B(c2, std::vector<float>(c1,0.0)); // c2 * c1

	init(A, c1, c2, key1);
	init(B, c2, c1, key2);
	auto r = mm(A,B, c1,c2);
	if(debug && key1 == 99)
		printf("compute method %f\n", r);
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
	std::map<long, bool> mapDisabledKeys;
	if(!myId){
		for(long i=0;i<nkeys; ++i) {
			map[i] = 0;
			mapDisabledKeys[i] = false;
		}
	}
	
	std::vector<float> V(nkeys, 0);
	std::vector<keypair> queueKeys;
	// bool resetkey1=false;
	// bool resetkey2=false;

	//KEYS DATATYPE
	int blockLengths_k[1] = {2};
	MPI_Aint lb_k, extent_k;
	MPI_Type_get_extent(MPI_LONG, &lb_k, &extent_k); 
	MPI_Aint disp_k[2] = {0, 2*extent_k};
	MPI_Datatype types_k[1] = {MPI_LONG}; 

	// Create the datatype for the parameters
	MPI_Datatype keysType;
	MPI_Type_create_struct(1, blockLengths_k, disp_k, types_k, &keysType);
	MPI_Type_commit(&keysType);

	//COMP DATATYPE
	int blockLengths_c[2] = {4,1};
	MPI_Aint lb_c, extent_c;
	MPI_Type_get_extent(MPI_LONG, &lb_c, &extent_c); 
	MPI_Aint disp_c[2] = {0, 4*extent_c};
	MPI_Datatype types_c[2] = {MPI_LONG, MPI_FLOAT}; 

	// Create the datatype for the parameters
	MPI_Datatype compType;
	MPI_Type_create_struct(2, blockLengths_c, disp_c, types_c, &compType);
	MPI_Type_commit(&compType);

	computation c;

	//MPI_Request req;
	MPI_Request rq_send;// rq_recv;
	//MPI_Status status;

	double start = MPI_Wtime();

	//START MAP FILLING
	// first section
	if(myId) { // i am a working node1

		bool check = true;
		while(check){
			// if(debug){
			// 	std::cout << "hello, i am the node " << myId << " and i am waiting to receive a message" << std::endl;
			// }
			//MPI_Irecv(&c, 1, compType, 0, COMPUTE_DATA, MPI_COMM_WORLD, &rq_recv);
			MPI_Recv(&c, 1, compType, 0, KEY_GENERATOR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if(debug && (c.m_1 != 0 || c.m_2 != 0)){
				printf("hello, i am the node %d and i have received the values: %ld, %ld, %ld, %ld\n", myId, c.key1, c.key2, c.m_1, c.m_2);
			}
			
			if(c.key1 == -1) {
				check = false;
				if(debug){
					std::cout << myId << " is leaving" << std::endl;
				}
			} else if (c.m_1 != 0 && c.m_2 != 0) {

				auto r = compute(c.m_1, c.m_2, c.key1, c.key2);
				V[c.key1] += r;
				c.result = r;
				c.key1 = c.key1;
				//MPI_Send(&k, 1, keysType, 0, COMPUTE_DATA, MPI_COMM_WORLD);
				if(c.key1 == 93){
					printf("hello, i am the node %d and i am computing : %f, %ld, %ld, %f \n", myId, c.result, c.key1, c.key2, r);
				}
				MPI_Isend(&c, 1, compType, 0, COMPUTE_DATA, MPI_COMM_WORLD, &rq_send);

				if(debug){
					std::cout << "V[key1] " << V[c.key1] << " with key: " << c.key1 << std::endl << "r: "<< r << std::endl;	
				}
			}	
		} 
		if(debug){
			for(long i=0;i<nkeys; ++i) {
				std::printf("key %ld : %f\n", i, V[i]);
			}
			std::cout << "END key generation for node " << myId << std::endl;
		} 
		
		// resetta 1
		// generazione key 1 -> lista k1 e k2 
	} else { // i am the server
		int roundrobin = 1;
		int disabledKeysNumber = 0;
		keypair tmpKey;
		keypair tmpKey2;
		int i = 0;
		while(i<length || disabledKeysNumber > 0 || !queueKeys.empty()) {
			if(roundrobin>(numP-1)){
				roundrobin = 1; 
				if(debug){
					std::cout << "roundrobin resettato a 1" << std::endl;
				}
			}

			// std::cout<< "disabledKeys sono più di 1" << std::endl;

			if(i<length){
				++i;
				key1 = random(0, nkeys-1);  // value in [0,nkeys[
				key2 = random(0, nkeys-1);  // value in [0,nkeys[
				if (key1 == key2) // only distinct values in the pair
					key1 = (key1+1) % nkeys;

				// qui devo controllare se una delle due chiavi è disabilitata se si metto in coda altrimenti aggiorno la mappa con tutte le conseguenze
				tmpKey.key1 = key1;
				tmpKey.key2 = key2;

				queueKeys.emplace_back(tmpKey);
				
				// printf("hello, i am the node server and i am generating: %ld, %ld last elem: %ld \n", key1, key2, (*(queueKeys.end())).key1);
			}
			if(disabledKeysNumber > 0) {
				MPI_Recv(&c, 1, compType, MPI_ANY_SOURCE, COMPUTE_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// MPI_Irecv(&k, 1, keysType, MPI_ANY_SOURCE, COMPUTE_DATA, MPI_COMM_WORLD, &rq_recv);
				if((c.key1 == 80 || c.key1 == 93) && false){
					printf("hello, i am the node server and i have received the values: %f, %ld\n", c.result, c.key1);
				}
				// c.m_1 //result
				// c.key1 //key to reset
				--disabledKeysNumber;
				auto _r1 = static_cast<unsigned long>(c.result) % SIZE;
				map[c.key1] = (_r1>(SIZE/2)) ? 0 : _r1;
			} else {
				if(!queueKeys.empty()) {
					tmpKey2 = *(queueKeys.begin());
					// qui controllo sempre che il primo non sia una coppia disabilitata.
					queueKeys.erase(queueKeys.begin());
					map[tmpKey2.key1]++;  // count the number of key1 keys
					map[tmpKey2.key2]++;  // count the number of key2 keys

					if (map[tmpKey2.key1] == SIZE && map[tmpKey2.key2] != 0) {
						c.m_1 = map[tmpKey2.key1];
						c.m_2 = map[tmpKey2.key2];
						c.key1 = tmpKey2.key1;
						c.key2 = tmpKey2.key2;
						++disabledKeysNumber;
						mapDisabledKeys[tmpKey2.key1] = true;
						if(debug){
							printf("hello, i am the node server and i have sent the values: %ld, %ld, %ld, %ld, round robin: %d \n", c.key1, c.key2, c.m_1, c.m_2, roundrobin);
						}
						MPI_Isend(&c, 1, compType, roundrobin, KEY_GENERATOR, MPI_COMM_WORLD, &rq_send);
						++roundrobin;
									
						if(roundrobin>(numP-1)){
							roundrobin = 1; 
							if(debug){
								std::cout << "roundrobin resettato a 1" << std::endl;
							}
						}
					}
					// if key2 reaches the SIZE limit ....
					if (map[tmpKey2.key2] == SIZE && map[tmpKey2.key1] != 0) {		
						c.m_1 = map[tmpKey2.key2];
						c.m_2 = map[tmpKey2.key1];
						c.key1 = tmpKey2.key2;
						c.key2 = tmpKey2.key1;
						++disabledKeysNumber;
						mapDisabledKeys[tmpKey2.key2] = true;
						MPI_Isend(&c, 1, compType, roundrobin, KEY_GENERATOR, MPI_COMM_WORLD, &rq_send);
						if(debug){
							printf("hello, i am the node server and i have sent the values: %ld, %ld, %ld, %ld, round robin: %d\n", c.key1, c.key2, c.m_1, c.m_2, roundrobin);
						}
						++roundrobin;
					}
				}
			} 
		} 
		c.key1 = -1;
		for(int g=1; g<numP; g++){
			if(debug){
				std::cout << "i am a working node " << myId << " and i am sending " << c.key1 << " to " << g << std::endl;
			}
			// MPI_Isend(&c, 1, compType, g, COMPUTE_DATA, MPI_COMM_WORLD, &rq_send);
			MPI_Send(&c, 1, compType, g, KEY_GENERATOR, MPI_COMM_WORLD);	
		} 
		if(debug){
			std::cout << "END key sharing and compute" << std::endl;
			for(long i=0;i<nkeys; ++i) {
				std::printf("key %ld : %f\n", i, V[i]);
			}

			for(long i=0;i<nkeys; ++i){
				std::printf("[MAP] key %ld : %ld\n", i, map[i]);
			}
		}
		// map[80] = 49;
		// map[93] = 10;
	}

	return 0;
	// second section
	// compute the last values in the other nodes but managements by server
	if(!myId) { // i am the server
		int roundrobin = 1; //id of the working node that we send the information for run the method compute
		if(debug2) {
			std::cout << "i am in the server" << std::endl;
		}
		for(long i=0;i<nkeys; ++i) {
			for(long j=0;j<nkeys; ++j) {
				if (i==j) continue;
				if (map[i]>0 && map[j]>0) {
					c.m_1 = map[i];
					c.m_2 = map[j];
					c.key1 = i;
					c.key2 = j;

					if(debug2 && c.key1 == 99){
						std::cout << "PRINT OF c1: key1 = " << c.key1 << " key2= " << c.key2 << " m1= " << c.m_1 << " m2= " << c.m_2 << std::endl;
					}

					// MPI_Isend(&c, 1, compType, roundrobin, COMPUTE_DATA, MPI_COMM_WORLD, &rq_send);
					MPI_Send(&c, 1, compType, roundrobin, COMPUTE_DATA, MPI_COMM_WORLD);
					++roundrobin;

					if(roundrobin>(numP-1)){
						roundrobin = 1; 
						if(debug2){
							std::cout << "roundrobin resettato a 1" << std::endl;
						}
					}
					
					c.m_1 = map[j];
					c.m_2 = map[i];
					c.key1 = j;
					c.key2 = i;

					if(debug2){
						std::cout << "PRINT OF c2: key1 = " << c.key1 << " key2= " << c.key2 << " m1= " << c.m_1 << " m2= " << c.m_2 << std::endl;
					}

					// MPI_Isend(&c, 1, compType, roundrobin, COMPUTE_DATA, MPI_COMM_WORLD, &rq_send);
					MPI_Send(&c, 1, compType, roundrobin, COMPUTE_DATA, MPI_COMM_WORLD);
					++roundrobin;

					if(roundrobin>(numP-1)){
						roundrobin = 1; 
						if(debug2 && false){
							std::cout << "roundrobin resettato a 1" << std::endl;
						}
					}
				}
			}
		} 
		// we have to send to all the working node a message of EOS for the computing phase
		c.key1 = -1;
		for(int g=1; g<numP; g++){
			if(debug2){
				std::cout << "i am a working node " << myId << " and i am sending " << c.key1 << " to " << g << std::endl;
			}
			// MPI_Isend(&c, 1, compType, g, COMPUTE_DATA, MPI_COMM_WORLD, &rq_send);
			MPI_Send(&c, 1, compType, g, COMPUTE_DATA, MPI_COMM_WORLD);	
		}
		
		std::vector<float> tempV(nkeys, 0);
		std::vector<std::vector<float>> tempVV;
		for(int j=0; j<(numP-1); j++){
			MPI_Recv(tempV.data(), nkeys, MPI_FLOAT, MPI_ANY_SOURCE, SEND_VECTOR_V, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//MPI_Irecv(tempV.data(), nkeys, MPI_FLOAT, MPI_ANY_SOURCE, SEND_VECTOR_V, MPI_COMM_WORLD, &rq_recv);
			tempVV.emplace_back(tempV);
			//MPI_Wait(&rq_recv, &status);
		}

		for(int j=0; j<(numP-1); j++) { // i have to merge the vector in only one
			// MPI_Irecv(tempV.data(), nkeys, MPI_FLOAT, MPI_ANY_SOURCE, SEND_VECTOR_V, MPI_COMM_WORLD, &req);
			// MPI_Recv(tempV.data(), nkeys, MPI_FLOAT, MPI_ANY_SOURCE, SEND_VECTOR_V, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for(int i=0; i<nkeys;i++){
				if(debug2){
					std::cout << "V[i]: " << V[i] <<" temp: "<<tempV[i] << std::endl;
				}
				if(!tempVV[j].empty()){
					V[i] += (tempVV[j])[i];
				}
				//tempV[i];
				if(debug2){
					std::cout << "V[i]: " << V[i] << std::endl;
				}
			}
		}
		
	} else { //we leave the server and enter in the nodes
		if(debug2){
			std::cout << "hello, i am the node " << myId << std::endl;
			std::cout << "and V[99] " << V[99] << std::endl;
		}
		bool check = true;
		while(check){
			if(debug2){
				std::cout << "hello, i am the node " << myId << " and i am waiting to receive a message" << std::endl;
			}
			//MPI_Irecv(&c, 1, compType, 0, COMPUTE_DATA, MPI_COMM_WORLD, &rq_recv);
			MPI_Recv(&c, 1, compType, 0, COMPUTE_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//MPI_Wait(&rq_recv, &status);
			if(debug2 && c.key1 == 99){
				printf("hello, i am the node %d and i have received the values: %ld, %ld, %ld, %ld\n", myId, c.key1, c.key2, c.m_1, c.m_2);
			}
			if(c.key1 == -1){
				check = false;
				if(debug2){
					std::cout << myId << " is leaving" << std::endl;
				}
			} else {
				auto r = compute(c.m_1, c.m_2, c.key1, c.key2);
				V[c.key1] += r;
				if(debug2 && c.key1 == 99){
					// std::cout << "V[key1] " << V[c.key1] << " with key: " << c.key1 << "\t r: "<< r << std::endl;	
					std::cout << "r: "<< r << std::endl;	
				}
			}	
		}
		if(debug2){
			std::cout << "hello, i am the node " << myId << " and i am sending to V" << std::endl;
		}
		// MPI_Isend(V.data(), V.size(), MPI_FLOAT, 0, SEND_VECTOR_V, MPI_COMM_WORLD, &rq_send);
		MPI_Send(V.data(), V.size(), MPI_FLOAT, 0, SEND_VECTOR_V, MPI_COMM_WORLD);
		
		if(debug2){
			for(int h=0; h<nkeys; h++){
				std::cout << h << " V[h]: " << V[h] << std::endl;
			}
		}
	}
	

	double end = MPI_Wtime();

	// printing the results
	if (print && !myId) {
		for(long i=0;i<nkeys; ++i)
			std::printf("key %ld : %f\n", i, V[i]);
	}
	std::cout << "Time with " << numP << " processes: " << end-start << " seconds" << std::endl;
	
	MPI_Type_free(&keysType);
	MPI_Type_free(&compType);

	// Terminate MPI
	MPI_Finalize();
	return 0;

}

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