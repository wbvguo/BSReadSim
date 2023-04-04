#include <cstdlib>
#include <random>
#include <chrono>
#include <iostream>


std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);
std::bernoulli_distribution dis_ber(0.5);

// generate a random number between 0 and 1
double random_num1() {
    return drand48();
}

double random_num2() {
    return dis(gen);
}

double random_num3() {
    return drand48()>0.5?0:1;
}

double random_num4() {
    return dis_ber(gen);
}



int main() {
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100000000; i++) { random_num1();}
    // stop the timer and calculate the duration
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // output the duration in microseconds
    std::cout << "RN1 Duration: " << duration.count() << " microseconds" << std::endl;


    start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100000000; i++) { random_num2();}
    // stop the timer and calculate the duration
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // output the duration in microseconds
    std::cout << "RN2 Duration: " << duration.count() << " microseconds" << std::endl;


    start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100000000; i++) { random_num3();}
    // stop the timer and calculate the duration
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // output the duration in microseconds
    std::cout << "RN3 Duration: " << duration.count() << " microseconds" << std::endl;


    start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100000000; i++) { random_num4();}
    // stop the timer and calculate the duration
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // output the duration in microseconds
    std::cout << "RN4 Duration: " << duration.count() << " microseconds" << std::endl;

    return 0;
}