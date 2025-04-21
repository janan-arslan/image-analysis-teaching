#include "RandomStrategy.hpp"
#include <algorithm>
#include <vector>
#include <chrono>
#include <random>
#include <iostream>

namespace sevens {

// Constructor seeds the RNG
RandomStrategy::RandomStrategy() {
    auto seed = static_cast<unsigned long>(
        std::chrono::system_clock::now().time_since_epoch().count()
    );
    rng.seed(seed);
}

void RandomStrategy::initialize(uint64_t playerID) {
    myID = playerID;
    // No special initialization needed for this simplistic version
}

int RandomStrategy::selectCardToPlay(
    const std::vector<Card>& hand,
    const std::unordered_map<uint64_t, std::unordered_map<uint64_t, bool>>& tableLayout)
{
    // Very simplified logic:
    // 1. If our hand is empty, we can't play => return -1
    // 2. Otherwise, pick a random index in the hand (no adjacency checks)

    if (hand.empty()) {
        return -1; // pass
    }

    // Uniform random index from 0 to hand.size()-1
    std::uniform_int_distribution<int> dist(0, static_cast<int>(hand.size()) - 1);
    int idx = dist(rng);
    return idx;
}

void RandomStrategy::observeMove(uint64_t /*playerID*/, const Card& /*playedCard*/) {
    // This simplified strategy ignores other players' moves
}

void RandomStrategy::observePass(uint64_t /*playerID*/) {
    // This simplified strategy ignores passes
}

std::string RandomStrategy::getName() const {
    return "RandomStrategy";
}

} // namespace sevens
