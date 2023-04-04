#include <cstdlib>
#include <random>
#include <chrono>
#include <iostream>

void rand_gen(std::mt19937 &gen){
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for(int i = 0; i < 10; i++){
        printf("%f\n", dis(gen));
    }
}


int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    printf("random\n");
    rand_gen(gen);

    
    printf("seed\n");
    std::mt19937 gen2 {1000};
    rand_gen(gen2);

    printf("seed2\n");
    rand_gen(gen2);
    return 0;
}