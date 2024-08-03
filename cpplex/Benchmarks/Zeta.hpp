#include <Complex.hpp>
#include <Constants.hpp>
#include <iostream>
#include <cassert>
#include <cmath>

using namespace cpplex;

void benchCpplexZetaSum(){
    const int N = 1e+9;
    Complex s = 0.5 + 14.1347251417346937904572519835625_j;
    Complex zetaTerm = 1/(1 - pow(2, 1 - s));

    Complex sum = 0; 
    for(int n = 1; n <= N; n++){
        Complex etaTerm = pow(-1, n-1) * 1 / (pow(n, s));

        sum += etaTerm * zetaTerm;
    }
    std::cout << sum << "\n";
}

void benchStdZetaSum(){
    const int N = 1e+9;
    std::complex<double> s = 0.5 + 14.1347251417346937904572519835625i;
    std::complex<double> zetaTerm = std::complex<double>(1)/(std::complex<double>(1) - pow(2, std::complex<double>(1) - s));


    std::complex<double> sum = 0; 
    for(int n = 1; n <= N; n++){
        std::complex<double> etaTerm = std::pow(-1, n-1) * 1 / (pow(n, s));

        sum += etaTerm * zetaTerm;
    }
    std::cout << sum << "\n";
}