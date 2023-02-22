#include <stdio.h>
#include <random>
#include <vector>

int MAX_SIZE = 1000000;
std::random_device rd;    
std::mt19937 gen(rd());

std::vector<float> efficiency_vec;

int main() {
    char *fname = "/home/wbguo/iproject/BSReadSim/test/data/weights_descend.txt";
    FILE* fp = fopen(fname, "r");
    if(fp==NULL){fprintf(stderr, "[%s] ERROR: open capture efficiency file: %s failed. Exit...", __func__, fname); exit (EXIT_FAILURE);}
    
    float eff_prob;

    while(fscanf(fp, "%f", &eff_prob) == 1) {
        // printf("%f\n", eff_prob);
        efficiency_vec.push_back(eff_prob);
    }

    std::discrete_distribution<int> dis(efficiency_vec.begin(), efficiency_vec.end());
    
    printf("%ld\n", sizeof(dis));

    for(int i =0; i< MAX_SIZE; ++i){
        //dis(gen);
        printf("%d\n", dis(gen));
    }
    return 0;
}

//g++