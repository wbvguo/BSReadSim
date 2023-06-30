#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

int main() {
  for(int i =0; i < 4; ++i){
    printf("i=%d\n", i);
    
    int offset[2] = {0, 0};
    printf("%d\n", offset[0]);
    offset[0] +=1;
    printf("%d\n", offset[0]);
  }
}
