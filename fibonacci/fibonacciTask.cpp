#include<iostream>
#include<thread>
#include<cstdint>
#include<future>
#include<vector>
#include<functional>

template<
    typename Func,
    typename ... Args,
    typename Rtrn=typename std::result_of<Func(Args...)>::type
    >
    auto make_task(
        Func && func,
        Args && ...args
    ) -> std::packaged_task<Rtrn(void)> {

        auto aux = std::bind(
            std::forward<Func>(func),
            std::forward<Args>(args)...
        );

        auto task = std::packaged_task<Rtrn(void)>(aux);
        return task;        
}

uint64_t fibo(uint64_t n) {
    uint64_t a_0 = 0, a_1 = 1;
    for(uint64_t index=0;index<n;index++) {
        const uint64_t tmp = a_0; a_0 = a_1; a_1 += tmp;   
    }
    return a_1;
}

int main(int argc, char const *argv[])
{
    const uint64_t num_threads = 32;
    
    std::vector<std::future<uint64_t>> results;
    std::vector<std::thread> threads;

    for(uint64_t index=0;index<num_threads;index++) {
        auto task = make_task(fibo, index);
        results.emplace_back(task.get_future());
        threads.emplace_back(std::move(task));
    }

    for(auto& result:results) std::cout<<result.get()<<std::endl;
    for(auto& thread:threads) thread.join();
}
