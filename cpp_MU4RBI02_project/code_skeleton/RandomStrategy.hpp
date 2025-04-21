#pragma once

#include "PlayerStrategy.hpp"
#include <random>
#include <chrono>

namespace sevens {

/**
 * A simple strategy that selects a random playable card.
 */
class RandomStrategy : public PlayerStrategy {
public:
    RandomStrategy();
    ~RandomStrategy() override = default;
    
    // PlayerStrategy interface
    void initialize(uint64_t playerID) override;
    int selectCardToPlay(
        const std::vector<Card>& hand,
        const std::unordered_map<uint64_t, std::unordered_map<uint64_t, bool>>& tableLayout) override;
    void observeMove(uint64_t playerID, const Card& playedCard) override;
    void observePass(uint64_t playerID) override;
    std::string getName() const override;
    
private:
    uint64_t myID;
    std::mt19937 rng;
};

} // namespace sevens
