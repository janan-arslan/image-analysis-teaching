#include <iostream>
#include <string>

// Include your framework files here...
// e.g. #include "MyGameMapper.hpp"
// #include "RandomStrategy.hpp"
// #include "GreedyStrategy.hpp"
// #include "StrategyLoader.hpp"

int main(int argc, char* argv[]) {
    // This is a minimal skeleton for demonstration purposes.
    // Students should integrate their classes or call the relevant
    // game logic from MyGameMapper (or other classes) as needed.
    
    if (argc < 2) {
        std::cout << "Usage: ./sevens_game [mode] [optional libs...]\n";
        return 1;
    }
    
    std::string mode = argv[1];
    
    if (mode == "internal") {
        std::cout << "[main] Internal mode is not fully implemented.\n";
        // TODO: Possibly create MyGameMapper instance and simulate 
        // with default random approach
    }
    else if (mode == "demo") {
        std::cout << "[main] Demo mode is not fully implemented.\n";
        // TODO: Possibly register built-in strategies and run a game
    }
    else if (mode == "competition") {
        std::cout << "[main] Competition mode is not fully implemented.\n";
        // TODO: Dynamically load .so libraries from argv[2..] 
        // and simulate a competition
    }
    else {
        std::cerr << "[main] Unknown mode: " << mode << std::endl;
    }
    
    return 0;
}
