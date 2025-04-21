#pragma once

#include "PlayerStrategy.hpp"
#include <memory>
#include <string>
#include <dlfcn.h>
#include <stdexcept>
#include <iostream>

namespace sevens {

/**
 * Utility class for loading player strategies from shared libraries.
 */
class StrategyLoader {
public:
    static std::shared_ptr<PlayerStrategy> loadFromLibrary(const std::string& libraryPath) {
        // TODO: dynamic loading, or stub
        // For a skeleton, you might just throw an exception:
        throw std::runtime_error("StrategyLoader::loadFromLibrary is not implemented in skeleton.");
    }
};

} // namespace sevens
