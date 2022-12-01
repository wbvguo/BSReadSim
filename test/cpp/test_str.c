#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int main(){
    char *str = "C,CGG;T,CGA";
    int i = 0;    
    while (str[i])
    {
        printf("%c\n", str[i]);
        ++i;
    }
    return 0;
}