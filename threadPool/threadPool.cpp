#include <threadPool.hpp>
#include <iostream>
#include <future>
#include <vector>

using namespace std;
auto square(const uint64_t x) {
    return x*x;
}

int main(int argc, char const *argv[])
{
    ThreadPool TP(8);

    const uint64_t num_task = 128;
    vector<future<uint64_t>> futures;

    for(uint64_t x=0;x<num_task;x++) {
        futures.emplace_back(TP.enqueue(square, x));
    }

    for(auto& future : futures) {
        cout << future.get() << endl;
    }
}
