#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    char c[] = "ls";
    char *p = c;
    int x = ((int)(strcmp(p,"+")!=0)<<1) | (int)(strcmp(p,"-")!=0);
    fprintf(stdout, "hello world\n");
    fprintf(stdout, "%d\n", x);
}

