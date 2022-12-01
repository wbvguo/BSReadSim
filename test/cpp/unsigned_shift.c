# include <stdio.h>

int main(){
    unsigned short x = 1;
    unsigned short y = x << 15;
    printf("%i\n", y);
    printf("%i\n", y >> 15);
    return 0;
}