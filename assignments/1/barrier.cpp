#include <thread>
#include <iostream>
#include <barrier>
#include <vector>

using namespace std;

barrier myBarrier{4};

void work(std::chrono::seconds w) {
	auto end = std::chrono::steady_clock::now() + w;
    while(std::chrono::steady_clock::now() < end);	 
}

void helloWorld(int id, int id2) {
    printf("<%d,%d> hello world\n", id, id2);
    work(std::chrono::seconds(5));
    myBarrier.arrive_and_wait();
}

int main(int argc, char const *argv[])
{
    vector<thread> threads;
    int numThread = 3;
    int numTask = 6;

    for(int index=0;index<numTask;index++) {
        for(int k=0;k<numThread;k++) {
            threads.emplace_back(
                helloWorld, k, index
            );
        }
        myBarrier.arrive_and_wait();
    }

    for(auto& thread : threads) thread.join();
}
