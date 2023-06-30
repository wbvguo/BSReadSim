#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

int main() {
  int i = -256;
  int16_t j = --i;
  printf("%d %d\n", i, j);
  return 0;
}
