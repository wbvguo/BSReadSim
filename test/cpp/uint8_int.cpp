#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <chrono>

int main(){
    int max_length = 1000000;
    int *tmp_seq[2];    	// sequence
    tmp_seq[1] = (int*)calloc(max_length, 4);
    tmp_seq[2] = (int*)calloc(max_length, 4);
    auto start_time = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < max_length; ++i)
    {
        tmp_seq[1][i] = 0;
        tmp_seq[2][i] = 0;
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    fprintf(stdout, "%ld us\n", duration.count());
    
    uint8_t *tmp_seq2[2];    	// sequence
    tmp_seq2[1] = (uint8_t*)calloc(max_length, 1);
    tmp_seq2[2] = (uint8_t*)calloc(max_length, 1);
    start_time = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < max_length; ++i)
    {
        tmp_seq2[1][i] = 0;
        tmp_seq2[2][i] = 0;
    }
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    fprintf(stdout, "%ld us\n", duration.count());
    
    return 0;
}