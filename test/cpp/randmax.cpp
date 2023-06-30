#include <cstdlib>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <ctime>  

int main() {
    srand(time(0)&0x7fffffff);

    int x, t;
    int count, total;
    count = 0;
    total = 1000000;
    for (size_t i = 0; i < total; i++)
    {
        t = rand();
        x = t > (RAND_MAX>>1);
        count += x;
    }
    printf("%d %d %d\n", total, count, RAND_MAX>>1);
    return 0;
}

