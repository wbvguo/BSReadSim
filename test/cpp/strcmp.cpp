#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>


int main() {
    char a[] = "";
    char b[] = "None";
    char c[] = "abc";
    bool x0 = strcmp(a, "None");
    bool x1 = strlen(a);
    bool x2 = strcmp(a, "None") && strlen(a);
    bool y0 = strcmp(b, "None");
    bool y1 = strlen(b);
    bool y2 = strcmp(b, "None") && strlen(b);
    bool z0 = strcmp(c, "None");
    bool z1 = strlen(c);
    bool z2 = strcmp(c, "None") && strlen(c);
    
    fprintf(stderr, "%d %d %d\n", x0, x1, x2);
    fprintf(stderr, "%d %d %d\n", y0, y1, y2);
    fprintf(stderr, "%d %d %d\n", z0, z1, z2);
}