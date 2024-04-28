#include <iostream>
#include <vector>

using namespace std;


void blockCyclicDataDistribution(const uint64_t& id, 
    const uint64_t num_threads, const uint64_t chunk_size,
    const uint64_t num_tasks, vector<uint64_t> source) {
	
    const uint64_t stride = num_threads*chunk_size;
	const uint64_t offset = id*chunk_size;

	for(uint64_t lower = offset; lower < num_tasks; lower += stride) {
		const uint64_t upper = std::min(lower+chunk_size, num_tasks);
        cout<<"lower: "<<lower<<" upper:"<<upper<<endl;

		for(uint64_t row = lower; row < upper; row++) {
			cout<<"row: "<<row<<endl;
		}
	}
}

void printMatrix(vector<int> &M, int N) {
    for(uint64_t k = 0; k< N; ++k) {  
			for(uint64_t i = 0; i< (N-k); ++i) { 
                cout << "Matrix [" << k << "," << i << "] " << M[i*N+(i+k)]<<endl;
            }
    }
}

int main(int argc, char const *argv[])
{
    // uint64_t sizeChunck = 3;
    // uint64_t numThreads = 3;
    
    // vector<uint64_t> source{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    
    // for(uint64_t index=0;index<numThreads;index++){
    //     blockCyclicDataDistribution(index, numThreads, sizeChunck, source.size(), source);
    //     cout<<"#################"<<endl;
    // }
    uint64_t N = 3;
    
	std::vector<int> M(N*N, -1);
    printMatrix (M, N);
    for(uint64_t k = 0; k< N; ++k) {  
			for(uint64_t i = 0; i< (N-k); ++i) { 
                M[i*N+(i+k)] = i;
            }
    }
    cout << "@@@@@@@@@@@@@@@" << endl;
    printMatrix (M, N);
}