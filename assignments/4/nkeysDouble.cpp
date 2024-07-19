#include <cstdio>
#include <random>
#include <map>
#include <vector>
#include <string>

const long   SIZE = 64;

long random(const int &min, const int &max) {
	static std::mt19937 generator(117);
	std::uniform_int_distribution<long> distribution(min,max);
	return distribution(generator);
};		


void init(auto& M, const long c1, const long c2, const long key) {
	for(long i=0;i<c1;++i)
		for(long j=0;j<c2;++j)
			M[i][j] = (key-i-j)/static_cast<double>(SIZE); //static_cast<float>(SIZE)
}

// matrix multiplication:  C = A x B  A[c1][c2] B[c2][c1] C[c1][c1]
// mm returns the sum of the elements of the C matrix
auto mm(const auto& A, const auto& B, const long c1,const long c2) {

	double sum{0}; //float
	
    for (long i = 0; i < c1; i++) {
        for (long j = 0; j < c1; j++) {
            auto accum = double(0.0); //float
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
double compute(const long c1, const long c2, long key1, long key2) { //float

	std::vector<std::vector<double>> A(c1, std::vector<double>(c2,0.0)); // c1 * c2 //std::vector<std::vector<float>> std::vector<float>(c2,0.0)
	std::vector<std::vector<double>> B(c2, std::vector<double>(c1,0.0)); // c2 * c1 std::vector<std::vector<float>>

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
	
	std::vector<double> V(nkeys, 0); //std::vector<float>
	std::vector<double> V2(nkeys, 0); //std::vector<float>
	bool resetkey1=false;
	bool resetkey2=false;
	int numReset = 0;
	for(int i=0;i<length; ++i) {
		key1 = random(0, nkeys-1);  // value in [0,nkeys[
		key2 = random(0, nkeys-1);  // value in [0,nkeys[
		
		if (key1 == key2) // only distinct values in the pair
			key1 = (key1+1) % nkeys; 

		map[key1]++;  // count the number of key1 keys
		map[key2]++;  // count the number of key2 keys

		double r1; //float
		double r2; //float
		// if key1 reaches the SIZE limit, then do the computation and then
		// reset the counter ....
		if (map[key1] == SIZE && map[key2]!=0) {			
			r1= compute(map[key1], map[key2], key1, key2);

			V[key1] += r1;  // sum the partial values for key1
			// if(key1==0){
			// 	printf("V[key1]:%f \n", V[key1]);
			// }
			resetkey1=true;			
		}
		// if key2 reaches the SIZE limit ....
		if (map[key2] == SIZE && map[key1]!=0) {	
			r2= compute(map[key2], map[key1], key2, key1);
			// if(key2==0){
			// 		printf("R2 compute results: m1=%ld, m2=%ld, key1=%ld, key2=%ld, r=%f \n", map[key1], map[key2], key1, key2, r2);
			// 	}
			V[key2] += r2;  // sum the partial values for key1
			// if(key2==0){
			// 		printf("V[key2]:%f \n", V[key2]);
			// 	}
			resetkey2=true;
		}
		if (resetkey1) {
			if(key1 == 0) {
				++numReset;
			}
			// updating the map[key1] initial value before restarting
			// the computation
			auto _r1 = static_cast<unsigned long>(r1) % SIZE;
			// if(key1 == 0) {
			// 	printf("hello, the computed value is %f key: %ld, %ld, numReset: %d\n", r1, key1, key2, numReset);
			// }	
			map[key1] = (_r1>(SIZE/2)) ? 0 : _r1;
			resetkey1 = false;
		}
		if (resetkey2) {
			if(key2 == 0) {
				++numReset;
			}
			// updating the map[key2] initial value before restarting
			// the computation
			auto _r2 = static_cast<unsigned long>(r2) % SIZE;
			// if(key2 == 0) {
			// 	printf("hello, the computed value is %f key: %ld, %ld, numReset: %d\n", r2, key2, key1, numReset);
			// }	
			map[key2] = (_r2>(SIZE/2)) ? 0 : _r2;
			resetkey2 = false;
		}
	}
	for(int h=0; h<nkeys; h++){
			printf("h: %d , V[h]:%f \n", h,V[h]);
		}
	// for(long i=0;i<nkeys; ++i){
	// 		std::printf("[MAP] key %ld : %ld\n", i, map[i]);
	// }
	
	// for(long i=0;i<nkeys; ++i){
	// 		std::printf("key %ld : %f\n", i, V[i]);
	// }
	// // return 0;
	

	
// return 0; and V[99] 1507517.875000 
// return 0; and V[99] 1.50752e+06
	// std::printf("and V[99] %f \n", V[99]);
	// compute the last values
	for(long i=0;i<nkeys; ++i) {
		for(long j=0;j<nkeys; ++j) {
			if (i==j) continue;
			if (map[i]>0 && map[j]>0) {
				auto r1= compute(map[i], map[j], i, j);
				
				printf("R1 compute results: m1=%ld, m2=%ld, key1=%ld, key2=%ld, r=%f \n", map[i], map[j], i, j, r1);
				
				auto r2= compute(map[j], map[i], j, i);
				printf("R2 compute results: m1=%ld, m2=%ld, key1=%ld, key2=%ld, r=%f \n", map[j], map[i], j, i, r2);
				
				// if(i == 98) {
				// 	printf("hello, the computed value is %f key: %ld, %ld\n", r1, i, j);
				// }	
				// if(j == 98) {
				// 	printf("hello, the computed value is %f key: %ld, %ld\n", r2, j, i);
				// }	
				// if(i == 98){
				// 	//printf("V[key1] %f, %ld, %ld, %ld, %ld\n", V[i], i, j, map[i], map[j]);
				// 	printf("r %f\n", r1);
				// } else if(j == 98){
				// 	printf("r %f\n", r2);
				// }
				V2[i] += r1;
				V2[j] += r2;
			}
		}
	}
	
	for(long i=0;i<nkeys; ++i) {
		printf("NODE V[h]: h: %ld , V[h]:%f \n", i, V2[i]);
	}

	for(long i=0;i<nkeys; ++i) {
		
		V[i] += V2[i];
	}
	for(long i=0;i<nkeys; ++i) {
		printf("NODE V[h]: h: %ld , V[h]:%f \n", i, V[i]);
	}
	

	// printing the results
	if (print) {
		for(long i=0;i<nkeys; ++i)
			std::printf("key %ld : %f\n", i, V[i]);
	}

	
}
