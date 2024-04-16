#include<iostream>
#include<cstdint>
#include<vector>
#include<future>

uint64_t fibo(uint64_t n) {
    uint64_t a_0 = 0, a_1 = 1;
    for(uint64_t index=0;index<n;index++) {
        const uint64_t tmp = a_0; a_0 = a_1; a_1 += tmp; 
    }
    return a_1;
}

int main(int argc, char const *argv[])
{
    /* code */
    return 0;
}
