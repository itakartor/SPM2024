#include <iostream>
#include <future>
#include <thread>

using namespace std;

int main(int argc, char const *argv[])
{
    std::promise<void> promise;
    
    auto shared_future = promise.get_future().share();

    auto students = [&] (int myId) -> void {
        shared_future.get();
        cout << "<" << myId << "> Time to make coffe!" << endl;
    };

    thread my_thread0(students, 0);
    thread my_thread1(students, 1);

    this_thread::sleep_for(2s);

    promise.set_value();

    my_thread0.join();
    my_thread1.join();
}
