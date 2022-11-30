#include <stdio.h>

#include <random>

void myfunc(auto ptr, std::mt19937 *gen){
    int x = (*ptr)(*gen);
    printf("%d %d %d\n", x, ptr, gen);
}


int main(){
    std::random_device rd;    
    std::mt19937 gen(rd());

    std::uniform_int_distribution<int> dis(2, 10);

    printf("%d \n", dis);
    std::uniform_int_distribution<int> *ptr;
    printf("%d %d %d %d\n", sizeof(ptr), sizeof(&dis), sizeof(dis), sizeof(gen));
    ptr = &dis;
    auto ptr2= &dis;
    printf("%d \n", dis(gen));
    printf("%d \n", (*ptr)(gen));
    printf("%d \n", (*ptr2)(gen));
    printf("%d %d %d\n", ptr, &dis, &ptr2, &gen);
    printf("%d %d %d %d\n", sizeof(ptr), sizeof(&dis), sizeof(dis), sizeof(gen));

    myfunc(ptr2, &gen);

    return 0;
}