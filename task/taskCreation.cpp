#include <iostream>
#include <future>
#include <functional>
#include <thread>
#include <vector>
#include <cstdint>



template<
    typename Func,
    typename ... Args, 
    typename Rtrn = typename std::result_of<Func(Args...)>::type
    >
    auto make_task(
        Func && func,
        Args && ...args) -> std::packaged_task<Rtrn(void)> {
        
        //funziona ausialiare tipo wrapper
        auto aux = std::bind(
            std::forward<Func>(func),
            std::forward<Args>(args)...
        );

        auto task = std::packaged_task<Rtrn(void)>(aux);

        return task;
    }

bool comp(float value, int64_t threshold) {
    return value < threshold;
}

int main(int argc, char const *argv[])
{
    auto task = make_task(comp, 199.3, 0); // 0 result because is not true
    auto future = task.get_future();

    std::thread thread(std::move(task));
    thread.detach();

    std::cout << future.get() << std::endl;
}
