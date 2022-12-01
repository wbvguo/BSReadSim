#include <stdio.h>
#include <random>

std::random_device rd;    
std::mt19937 gen(rd());

void myfunc(auto ptr){
    int x = (*ptr)(gen);
    printf("%d %d %d %d\n", x, ptr, gen, &gen);
}

void myfunc2(auto ptr){
    int x = (*ptr)(gen);
    printf("%d %d %d %d\n", x, ptr, gen, &gen);
}

int main(){
    std::uniform_int_distribution<int> dis(2, 10);
    printf("%d %d\n", dis, gen);
    std::uniform_int_distribution<int> *ptr;
    printf("%d %d %d %d\n", sizeof(ptr), sizeof(&dis), sizeof(dis), sizeof(gen));
    ptr = &dis;
    auto ptr2= &dis;
    printf("%d \n", dis(gen));
    printf("%d \n", (*ptr)(gen));
    printf("%d \n", (*ptr2)(gen));
    printf("%d %d %d %d %d\n", ptr, &dis, ptr2, gen, &gen);
    printf("%d %d %d %d\n", sizeof(ptr), sizeof(&dis), sizeof(dis), sizeof(gen));

    myfunc(ptr);

    std::discrete_distribution<int> dis2({1, 2, 3, 4});
    auto ptr3 = &dis2;
    printf("%d %d\n", dis2, gen);
    printf("%d %d %d %d %d\n", ptr3, dis2, &dis2, gen, &gen);

    myfunc2(ptr3);
    return 0;
}