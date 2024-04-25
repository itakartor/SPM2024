#include <iostream>
#include <vector>
#include <deque>
#include <mutex>
#include <condition_variable>
#include <random>
#include <thread>
#include <stop_token>

using namespace std;

deque<int> dataq;

void producer(const stop_token *stoken, int id, condition_variable_any *cv, mutex *mtx) {
    int i = 0;
    while (!(*stoken).stop_requested()) {
        {
            lock_guard<mutex> lock(*mtx);
            dataq.push_back(i++);
            cout<<"Producer "<<id<< " has makes " << i <<endl;
        }
        (*cv).notify_one();
        this_thread::sleep_for(chrono::milliseconds(random()%100 + 1));
    }
    cout << "Producer " << id << " exits" << endl;
}

void consumer (const stop_token *stoken, int id, condition_variable_any *cv, mutex *mtx) {
    while(!(*stoken).stop_requested() || !dataq.empty()) {
        unique_lock<mutex> lock(*mtx);
        (*cv).wait(lock, *stoken, [] {return !dataq.empty();});
        if(!dataq.empty()) {
            uint64_t data = dataq.front();
            cout<<"Consumer "<<id<< " has consumed " << data <<endl;
            dataq.pop_front();
            lock.unlock();
        }
        this_thread::sleep_for(chrono::milliseconds(random()%100 + 1));
    }
    cout << "Consumer " << id << " exits" << endl;
}
int main(int argc, char const *argv[])
{
    int nConsumer = 2;
    int nProducer = 4;
    vector<jthread> producers;
    vector<jthread> consumers;

    stop_source stopSrc;
    stop_token stoken = stopSrc.get_token();
    condition_variable_any cv;

    mutex mutex;
    for(int i=0;i<nConsumer;i++) {
        consumers.emplace_back(consumer, &stoken, i, &cv, &mutex);
    }

    for(int i=0;i<nProducer;i++) {
        producers.emplace_back(producer, &stoken, i, &cv, &mutex);
    }

    this_thread::sleep_for(chrono::seconds(3));
    cout << "stopping main thread";
    stopSrc.request_stop();

    // not needed, here to make valgrind happy!
	for (auto& thread : consumers) thread.join();
	for (auto& thread : producers) thread.join();
}
