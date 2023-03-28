#include <stdio.h>

int main() {
    fprintf(stdout, "arr1");
    putc('\n', stdout);
    int ascii_hex[] =  {0x00, 0x10, 0x30, 0x70, 0x90, 0xb0, 0xf0, 
                        0x01, 0x11, 0x31, 0x71, 0x91, 0xb1, 0xf1, 
                        0x03, 0x13, 0x33, 0x73, 0x93, 0xb3, 0xf3};
    int ascii_hex2[]=  {0x00, 0x01, 0x03, 0x07, 0x09, 0x0b, 0x0f, 
                        0x10, 0x11, 0x13, 0x17, 0x19, 0x1b, 0x1f, 
                        0x30, 0x31, 0x33, 0x37, 0x39, 0x3b, 0x3f};
    int size = sizeof(ascii_hex) / sizeof(ascii_hex[0]);
    int ii=0;
    for(int i = 0; i < size; i++) {
        putc(ascii_hex[i], stdout);
        ++ii;
    }
    putc('\n', stdout);
    fprintf(stdout, "arr1 size %d", ii);
    putc('\n', stdout);
    
    fprintf(stdout, "arr2");
    putc('\n', stdout);
    ii=0;
    for(int i = 0; i < size; i++) {
        putc(ascii_hex2[i], stdout);
        ++ii;
    }
    putc('\n', stdout);
    fprintf(stdout, "arr2 size %d", ii);
    putc('\n', stdout);

    putc('\n', stdout);
    for (size_t i = 0; i < 10; i++)
    {
        fputc(i+33, stdout);
    }
    
    putc('\n', stdout);

}