#include "MyGameMapper.hpp"
#include <iostream>
#include <stdexcept>

namespace sevens {

MyGameMapper::MyGameMapper() {
    // TODO: Possibly seed random engine, etc.
}

void MyGameMapper::read_cards(const std::string& filename) {
    // TODO: e.g., call MyCardParser or implement logic here
    std::cout << "[MyGameMapper::read_cards] TODO: load or generate 52 cards.\n";
}

void MyGameMapper::read_game(const std::string& filename) {
    // TODO: e.g., initialize table_layout with 7s
    std::cout << "[MyGameMapper::read_game] TODO: set up table layout.\n";
}

bool MyGameMapper::hasRegisteredStrategies() const {
    // TODO: Return whether any strategies have been registered
    return false;
}

void MyGameMapper::registerStrategy(uint64_t playerID, std::shared_ptr<PlayerStrategy> strategy) {
    // TODO: store the strategy so we can use it during simulation
    (void)playerID;
    (void)strategy;
    std::cout << "[MyGameMapper::registerStrategy] TODO: store player's strategy.\n";
}

std::vector<std::pair<uint64_t, uint64_t>>
MyGameMapper::compute_game_progress(uint64_t numPlayers)
{
    // TODO: implement a quiet simulation
    // Return e.g. (playerID, finalRank) for each player
    std::cout << "[MyGameMapper::compute_game_progress] TODO.\n";
    return {};
}

std::vector<std::pair<uint64_t, uint64_t>>
MyGameMapper::compute_and_display_game(uint64_t numPlayers)
{
    // TODO: implement a verbose simulation
    std::cout << "[MyGameMapper::compute_and_display_game] TODO.\n";
    return {};
}

std::vector<std::pair<std::string, uint64_t>>
MyGameMapper::compute_game_progress(const std::vector<std::string>& playerNames)
{
    // Optional overload for name-based players
    (void)playerNames;
    std::cout << "[MyGameMapper::compute_game_progress(names)] TODO.\n";
    return {};
}

std::vector<std::pair<std::string, uint64_t>>
MyGameMapper::compute_and_display_game(const std::vector<std::string>& playerNames)
{
    // Optional overload for name-based players
    (void)playerNames;
    std::cout << "[MyGameMapper::compute_and_display_game(names)] TODO.\n";
    return {};
}

} // namespace sevens
