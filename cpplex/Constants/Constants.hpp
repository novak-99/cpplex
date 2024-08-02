#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "Complex/Complex.hpp"

namespace cpplex{
    const Complex INF = std::numeric_limits<double>::infinity();
    const Complex NINF = Complex(-std::numeric_limits<double>::infinity(), 0);

    const Complex li2 = 1.0451637801174927848445888891946131365226155781512015758329091440; 
    const Complex eulerGamma = 0.5772156649015328606065120900824024310421; 

}

#endif // CONSTANTS_HPP