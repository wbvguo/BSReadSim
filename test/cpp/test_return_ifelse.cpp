#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int sub(int k){
    return k == 0 ? 0 : 1;
}

int main(){
    int k;
    k = sub(0);
    fprintf(stderr, "k = %d\n", k);
    k = sub(100);
    fprintf(stderr, "k = %d\n", k);
    return 0;
}