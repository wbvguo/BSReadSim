#include <stdio.h>
#include <stdlib.h>

int main () {
   char str1[10], str2[10], str3[10];
   int year;
   FILE * fp;

   fp = fopen ("file.txt", "w+");
   fputs("We are in 2014 2 ", fp);
   
   rewind(fp);
   int x = fscanf(fp, "%s %s %s %d", str1, str2, str3, &year);
   
   printf("Read String1 |%s|\n", str1 );
   printf("Read String2 |%s|\n", str2 );
   printf("Read String3 |%s|\n", str3 );
   printf("Read Integer |%d|\n", year );
   printf("number of columns |%d|\n", x ); // output 4 or less, specified 4 variables

   fclose(fp);
   
   return(0);
}