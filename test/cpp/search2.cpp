#include <iostream>     // std::cout
#include <algorithm>    // std::search_n
#include <vector>       // std::vector

int main() {
    int a[] = { 1,2,3,4,4,4,1,2,3,4,4,4 };
    std::vector<int> find_vec = { 1, 2, 3 };

    int *ptr_start = std::begin(a);
    int *ptr_end = std::end(a);
    int *iter = std::begin(a);

    while (iter != ptr_end)
    {
        iter = std::search(ptr_start, ptr_end, find_vec.begin(), find_vec.end());
        std::cout << "one：" << iter - a << ",*iter = " << *iter << std::endl;
        ptr_start = iter + 3;
    }
    return 0;
}