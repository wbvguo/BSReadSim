# include <stdio.h>
# include <stdint.h>

int main(){
    uint64_t x = 2147483650;
    printf("%ld %d\n", x, (int)x);
    return 0;
}