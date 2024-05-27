#include <cstdio>
#include <random>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include "mpi.h"



int main(int argc, char const *argv[])
{
    MPI_Init(&argc, &argv);

	int myId;
	int numP;

    std::vector<float> V(4, 0);
    std::vector<float> localV(4, 0);
	// Get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &myId); 
    MPI_Comm_size(MPI_COMM_WORLD, &numP);

    if(myId) {
        for(int i=0;i<4;i++) {
            localV[i] = myId + i;
        }
        for(int i=0;i<4;i++) {
            std::cout<<"myId: "<<myId<<" localV: "<<localV[i]<<std::endl;
        }
    }

	MPI_Barrier(MPI_COMM_WORLD); //every process has A and B
    
    if(!myId) {
        // Gather results from all processes
    MPI_Gather(localV.data(),  // sendbuf 
			   4,
			   MPI_FLOAT,	
			   V.data(),	   // recvbuf 
			   4,         // recvcount; elements from each process
			   MPI_FLOAT,
			   0,              // target
			   MPI_COMM_WORLD);
        for(int i=0;i<4;i++) {
            std::cout<<"myId: 0"<<" V: "<<V[i]<<std::endl;
        }
    }
}
