#include <stdio.h>
#include <random>

std::random_device rd;    
std::mt19937 gen(rd());

int MAX_SIZE = 1000000000;

int main(int argc, const char** argv) {
    std::uniform_int_distribution<int> dis(0, MAX_SIZE);
    printf("%ld\n", sizeof(dis));
    for(int i =0; i< MAX_SIZE; ++i){
        dis(gen);
    }
    return 0;
}