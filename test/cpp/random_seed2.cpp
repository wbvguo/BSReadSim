#include <iostream>
#include <stdlib.h>
#include <time.h>

int main() {
    srand48(42);
    std::srand(time(0)&0x7fffffff);
        
    double rand_norm1 = drand48();
    int rand_norm2 = std::rand();
    fprintf(stdout, "rand_norm1: %f, rand_norm2: %d\n", rand_norm1, rand_norm2);

    return 0;
}
