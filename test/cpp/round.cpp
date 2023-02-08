#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <vector>

int main(){
    float a[] = {0,0.1,0.2,0.3,0.4,0.45,0.49,0.5,0.51,0.6,0.7,0.8,0.9,1.0};

    for (size_t i = 0; i < sizeof(a)/sizeof(a[0]); i++)
    {
        printf("%d\n", (int)(a[i]+0.5));
    }
}