#pragma once

#include "PlayerStrategy.hpp"

namespace sevens {

/**
 * A (placeholder) greedy strategy skeleton.
 */
class GreedyStrategy : public PlayerStrategy {
public:
    GreedyStrategy() = default;
    ~GreedyStrategy() override = default;
    
    void initialize(uint64_t playerID) override;
    int selectCardToPlay(
        const std::vector<Card>& hand,
        const std::unordered_map<uint64_t, std::unordered_map<uint64_t, bool>>& tableLayout) override;
    void observeMove(uint64_t playerID, const Card& playedCard) override;
    void observePass(uint64_t playerID) override;
    std::string getName() const override;
    
private:
    uint64_t myID;
};

} // namespace sevens
