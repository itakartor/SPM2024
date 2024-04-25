#include <iostream>
#include <threadPool.hpp>
#include <future>

using namespace std;
int main(int argc, char const *argv[])
{
    ThreadPool TP(4);

    auto square = [](const uint64_t x) {
        return x;
    };

    const uint64_t num_nodes = 13;
    vector<future<uint64_t>> futures;

    typedef function<void(uint64_t)> traverse_t;
    traverse_t traverse = [&] (uint64_t node) {
        if(node < num_nodes) {
            auto future = TP.enqueue(square, node);
            futures.emplace_back(move(future));

            traverse(2*node+1);
            traverse(2*node+2);
        }
    };

    traverse(0); //start from root;

    for(auto& future : futures) {
        cout << future.get() << endl;
    }
}
