#pragma once

#include "Generic_game_parser.hpp"
#include <string>
#include <vector>
#include <utility>

namespace sevens {

/**
 * Interface for simulating a Sevens game:
 *   - compute_game_progress(numPlayers)
 *   - compute_and_display_game(numPlayers)
 *
 * Subclasses must implement these pure virtual methods to run the logic.
 */
class Generic_game_mapper : public Generic_game_parser {
public:
    virtual std::vector<std::pair<uint64_t, uint64_t>>
    compute_game_progress(uint64_t numPlayers) = 0;

    virtual std::vector<std::pair<uint64_t, uint64_t>>
    compute_and_display_game(uint64_t numPlayers) = 0;

    // Overloads if you want name-based players
    virtual std::vector<std::pair<std::string, uint64_t>>
    compute_game_progress(const std::vector<std::string>& playerNames) {
        (void)playerNames; // avoid unused-parameter warning
        throw("Not implemented");
    }

    virtual std::vector<std::pair<std::string, uint64_t>>
    compute_and_display_game(const std::vector<std::string>& playerNames) {
        (void)playerNames;
        throw("Not implemented");
    }
};

} // namespace sevens
