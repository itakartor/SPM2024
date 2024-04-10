#include <cstdint> // uint64_t
#include <iostream> // std::cout std::endl
#include <vector> // std::vector
#include <thread> // std::thread


void say_hello_world(uint64_t idThread) {
    std::cout << "Hello from thread: " << idThread << std::endl;
}

int main(int argc, char const *argv[])
{
    const uint64_t num_thread{(argc == 2) ? std::stoul(argv[1]):10};
    std::vector<std::thread> threads;

    for (uint64_t id = 0; id < num_thread; id++) {
        // threads.push_back(
        //     std::thread(say_hello_world, id)
        // );

        // è meglio questo metodo per una questione di performance
        threads.emplace_back(
            say_hello_world,
            id
        );
    }

    for(auto& thread: threads) {
        thread.join();
    }
}

