#include <stdio.h>
#include <random>

std::random_device rd;    
std::mt19937 gen(rd());

int MAX_SIZE = 1000000;


int main(int argc, const char** argv) {
    std::discrete_distribution<int> dis(MAX_SIZE, 1, 1, std::bind2nd(std::plus<double>(),5.0));
    printf("%ld\n", sizeof(dis));
    for(int i =0; i< MAX_SIZE; ++i){
        printf("%d\n", dis(gen));
    }
    return 0;
}

//g++