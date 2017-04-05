//
// Created by phuschke on 3/31/17.
//

#pragma once

#include <ostream>
#include <iostream>
#include <vector>

template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    for (auto const& i : v)
        os << i << "\n";

    return os;
}