#include <stdio.h>
#include <stdlib.h>

int main() {
    char c[] = "+sdasd";
    char *p = c;
    fprintf(stdout, "hello world\n");
    fprintf(stdout, "%f\n", atof(p));
    return 0;
}