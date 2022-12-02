#include <stdio.h>
#include <random>
#include <vector>

std::random_device rd;    
std::mt19937 gen(rd());

class Solution{

};


int main(int argc, const char** argv) {
    std::discrete_distribution<int> *ptr1;
    std::uniform_int_distribution<int> *ptr2;
    std::discrete_distribution<int> dis({0,1,2,3,4,5});
    std::uniform_int_distribution<int> dis(0,5);

    return 0;
}