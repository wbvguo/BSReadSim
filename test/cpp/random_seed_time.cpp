#include <cstdlib>
#include <random>
#include <chrono>
#include <iostream>

void rand_gen(std::mt19937 &gen){
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for(int i = 0; i < 1000000; i++){
        float x = dis(gen);
    }
}

void rand_gen2(){
    for(int i = 0; i < 1000000; i++){
        float x = drand48();
    }
}

void rand_gen3(){
    for(int i = 0; i < 1000000; i++){
        float x = static_cast<float>(std::rand())/RAND_MAX;
    }
}


int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    auto start_time = std::chrono::high_resolution_clock::now();
    rand_gen(gen);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    fprintf(stdout, "rand3: %ld us\n", duration.count());

    start_time = std::chrono::high_resolution_clock::now();
    rand_gen2();
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    fprintf(stdout, "rand3: %ld us\n", duration.count());
    
    start_time = std::chrono::high_resolution_clock::now();
    rand_gen3();
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    fprintf(stdout, "rand3: %ld us\n", duration.count());

    fprintf(stdout, "rand3 value: %f %d\n", static_cast<float>(std::rand())/RAND_MAX, RAND_MAX);
}