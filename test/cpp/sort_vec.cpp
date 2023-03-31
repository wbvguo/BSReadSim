#include <algorithm> // for std::sort
#include <vector>
#include <stdio.h>

struct MyStruct {
    int value1;
    int value2;
};

int main() {
    std::vector<MyStruct> myVector = {
        {1, 2},
        {3, 4},
        {2, 3}
    };

    auto compareByValue1 = [](const MyStruct& a, const MyStruct& b) {
        return a.value1 < b.value1;
    };

    // Sort the vector based on value1
    std::sort(myVector.begin(), myVector.end(), compareByValue1);

    // Now the vector is sorted by value1
    for (size_t i = 0; i < myVector.size(); i++)
    {
        fprintf(stdout, "%d %d\n", myVector[i].value1, myVector[i].value2);
    }
}

