#include <cstdlib>
#include <random>
#include <chrono>
#include <iostream>

void rand_gen3(){
    for(int i = 0; i < 1000000; i++){
        float x = static_cast<float>(std::rand())/RAND_MAX;
    }
}


int main() {
    fprintf(stdout, "rand3 value: %f %d\n", static_cast<float>(std::rand())/RAND_MAX, RAND_MAX);
}