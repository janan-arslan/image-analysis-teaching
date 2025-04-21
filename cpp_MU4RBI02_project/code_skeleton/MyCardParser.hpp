#pragma once

#include "Generic_card_parser.hpp"

namespace sevens {

/**
 * Derived class from Generic_card_parser.
 * Subclasses must fill in read_cards(...) to create a deck or read from file.
 */
class MyCardParser : public Generic_card_parser {
public:
    MyCardParser() = default;
    ~MyCardParser() = default;

    void read_cards(const std::string& filename) override;
};

} // namespace sevens
