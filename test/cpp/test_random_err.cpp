#include <iostream>
#include <cstdlib>
#include <ctime>
#include <unordered_map>

int main() {
    // Seed the random number generator with the current time
    std::srand(static_cast<unsigned int>(std::time(0)));

    std::unordered_map<int, std::unordered_map<int, int>> valueCounts;

    // Repeat the process multiple times to get statistics
    int numIterations = 10000;

    for (int i = 0; i < numIterations; ++i) {
        for (int x = 0; x <= 3; ++x) {
            // Generate a random number between 0 and 2
            int randomMutation = std::rand() % 3;

            // Add 1 to the random number to ensure it's different from x
            int mutatedX = (x + randomMutation + 1)&3;

            valueCounts[x][mutatedX]++;
        }
    }

    // Print the value count statistics
    for (int x = 0; x <= 3; ++x) {
        std::cout << "Value " << x << " mutated to:" << std::endl;
        for (int mutatedX = 0; mutatedX <= 3; ++mutatedX) {
            int count = valueCounts[x][mutatedX];
            std::cout << "    " << mutatedX << ": " << count << " times" << std::endl;
        }
    }

    return 0;
}