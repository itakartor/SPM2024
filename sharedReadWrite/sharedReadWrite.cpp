#include <shared_mutex>
#include <iostream>
#include <thread>
#include <mutex>
#include <vector>

using namespace std;

void reader(int count, int id, int *shared_counter, shared_mutex *sMutex) {
    for(int i=0;i<count; ++i) {
        shared_lock<shared_mutex> lock(*sMutex);
        printf("Reader<%d> has read %d\n", id, *shared_counter);
    }
}

void writer(int count, int id, int *shared_counter, shared_mutex *sMutex) {
    for(int i=0;i<count; ++i) {
        unique_lock<shared_mutex> lock(*sMutex);
        ++(*shared_counter);
        printf("Writer<%d> has written %d\n", id, *shared_counter);
    }
}


int main(int argc, char const *argv[])
{
    int shared_counter = 0;
    shared_mutex sMutex;

    vector<thread> writers;
    for(int i=0;i<3; ++i) {
        writers.emplace_back(writer, 3, i, &shared_counter, &sMutex);
    }
    
    vector<thread> readers;
    for(int i=0;i<5; ++i) {
        readers.emplace_back(reader, 5, i, &shared_counter, &sMutex);
    }

    for (auto& thread : writers) thread.join();
	for (auto& thread : readers) thread.join();
}