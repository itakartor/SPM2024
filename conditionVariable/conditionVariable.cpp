#include <chrono>
#include <condition_variable>
#include <iostream>
#include <thread>

using namespace std::chrono_literals;

void student(int name, std::mutex *mutex, bool *time_for_breakfast, std::condition_variable *cv) {
    std::unique_lock<std::mutex> unique_lock(*mutex);

    while(!(*time_for_breakfast)) {
        if(!(*time_for_breakfast))std::cout <<"<"<< name << ">i have to wait" << std::endl;
        (*cv).wait(unique_lock);
    }
    std::cout <<"<"<< name << ">time to make coffe!" << std::endl;
    (*cv).notify_one();
}

int main(int argc, char const *argv[])
{
    std::mutex mutex;
    std::condition_variable cv;

    bool time_for_breakfast = false;

    // auto student = [&] () -> void {
    //     //this is for access to the lock
    //     std::unique_lock<std::mutex> unique_lock(mutex);

    //     while(!time_for_breakfast) {
    //         cv.wait(unique_lock);
    //     }
    //     std::cout << "time to make coffe!" << std::endl;

    // };
    std::thread my_thread2(student, 2, &mutex, &time_for_breakfast, &cv);
    std::this_thread::sleep_for(2s);
    std::thread my_thread1(student, 1, &mutex, &time_for_breakfast, &cv);
    std::this_thread::sleep_for(3s);
    

    {
        std::lock_guard<std::mutex> lock_guard(mutex);
        time_for_breakfast = true;
    }

    cv.notify_one();
    my_thread1.join();
    my_thread2.join();
}
