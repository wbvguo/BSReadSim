#include <iostream>
#include <map>

typedef struct frag_rrbs_rec {
    int pos_l, pos_r, start2;
    int8_t haplo = 0;
    int8_t strand = 0;
    int8_t cut_l = -1;
    int8_t cut_r = -1;
    std::map<int8_t, int8_t> n_cuts;

    bool operator<(const frag_rrbs_rec& other) const {
        if (pos_l != other.pos_l) {
            return pos_l < other.pos_l;
        } else {
            return pos_r < other.pos_r;
        }
    }
} frag_rrbs_rec;


typedef struct frag_rrbs_rec2 {
    int pos_l, pos_r, start2;
    int8_t haplo = 0;
    int8_t strand = 0;
    int8_t cut_l = -1;
    int8_t cut_r = -1;
    std::map<int8_t, int8_t> n_cuts;
} frag_rrbs_rec2;



int main() {
    // create a map with frag_rrbs_rec as key and int as value
    std::map<frag_rrbs_rec, int> frag_map;

    // add some data to the map
    frag_rrbs_rec frag1{10, 20, 30};
    frag_rrbs_rec frag2{15, 25, 35};
    frag_rrbs_rec frag3{20, 30, 40};
    frag_rrbs_rec frag4{10, 30, 50};

    frag_rrbs_rec2 frag5{10, 20, 30};

    frag_map[frag1] = 1;
    frag_map[frag2] = 2;
    frag_map[frag3] = 3;
    frag_map[frag4] = 4;

    // iterate through the map and print the keys and values
    for (auto& p : frag_map) {
        std::cout << "Key: (" << p.first.pos_l << "," << p.first.pos_r << "," << p.first.start2 << ") "
                  << "Value: " << p.second << std::endl;
    }

    fprintf(stdout, "frag4: %ld, frag5: %ld \n", sizeof(frag4), sizeof(frag5));

    return 0;
}