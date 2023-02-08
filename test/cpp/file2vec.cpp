#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <vector>

#define MAX_BIN 128

std::vector<float> efficiency_vec;

int main() {
    char *fname = "/home/wbguo/iproject/BSReadSim/test/data/eff.txt";
    FILE* fp = fopen(fname, "r");
    if(fp==NULL){fprintf(stderr, "[%s] ERROR: open capture efficiency file: %s failed. Exit...", __func__, fname); exit (EXIT_FAILURE);}
    
    float eff_prob;

    while(fscanf(fp, "%f", &eff_prob) == 1) {
        // printf("%f\n", eff_prob);
        efficiency_vec.push_back(eff_prob);
    }

    for (size_t i = 0; i < efficiency_vec.size(); i++)
    {
        printf("%ld %f\n", i, efficiency_vec[i]);
    }

    fclose(fp);
}
