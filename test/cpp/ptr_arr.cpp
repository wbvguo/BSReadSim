#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <vector>

int main(){
    int a[] = {0,0,0,0, 1,1,2,2, 3,2,1,0, 0,0,0,0, 3,1,2,0};

    int *p1 = a;
    int *p2 = *(&a+1);
    int *p3 = std::begin(a);
    int *p4 = std::end(a);

    printf("%d %d\n", p1,p2);
    printf("%d %d\n", p3,p4);

    printf("%d \n",  &a);
    printf("%d \n",  (&a+1));
}