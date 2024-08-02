#include <Complex/Complex.hpp>
#include <complex>

using namespace cpplex;
using namespace std::complex_literals; // for C++'s literals

void benchCpplexAddition(){
    const int N = 1e+9;

    Complex z1 = 1 + 1_j;
    Complex z2 = 3 + 3_j;
    Complex w = 0; 

    for(int i = 0; i < N; i++){
        w += z1 + z2; 
    }

    std::cout << w << "\n";
}

void benchStdAddition(){
    const int N = 1e+9;

    std::complex<double> z1 = 1.0 + 1i;
    std::complex<double> z2 = 3.0 + 3i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += z1 + z2; 
    }
    std::cout << w << "\n";
}

void benchCpplexSubtraction(){
    const int N = 1e+9;

    Complex z1 = 1 + 1_j;
    Complex z2 = 3 + 3_j;
    Complex w = 0; 

    for(int i = 0; i < N; i++){
        w += z1 - z2; 
    }

    std::cout << w << "\n";
}

void benchStdSubtraction(){
    const int N = 1e+9;

    std::complex<double> z1 = 1.0 + 1i;
    std::complex<double> z2 = 3.0 + 3i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += z1 - z2; 
    }
    std::cout << w << "\n";
}

void benchCpplexMultiplication(){
    const int N = 1e+9;

    Complex z1 = 1 + 1_j;
    Complex z2 = 3 + 3_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += z1 * z2; 
    }
    std::cout << w << "\n";
}

void benchStdMultiplication(){
    const int N = 1e+9;

    std::complex<double> z1 = 1.0 + 1i;
    std::complex<double> z2 = 3.0 + 3i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += z1 * z2; 
    }
    std::cout << w << "\n";
}

void benchCpplexDivision(){
    const int N = 1e+9;

    Complex z1 = 1 + 1_j;
    Complex z2 = 3 + 3_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += z1 / z2; 
    }
    std::cout << w << "\n";
}

void benchStdDivision(){
    const int N = 1e+9;

    std::complex<double> z1 = 1.0 + 1i;
    std::complex<double> z2 = 3.0 + 3i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += z1 / z2; 
    }
    std::cout << w << "\n";
}

void benchCpplexMod(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += mod(z);
    }
    std::cout << w << "\n";
}

void benchStdAbs(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::abs(z); 
    }
    std::cout << w << "\n";
}

void benchCpplexArg(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += arg(z);
    }
    std::cout << w << "\n";
}

void benchStdArg(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::arg(z); 
    }
    std::cout << w << "\n";
}

void benchCpplexConj(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += conj(z);
    }
    std::cout << w << "\n";
}

void benchStdConj(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::conj(z); 
    }
    std::cout << w << "\n";
}

void benchCpplexProj(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += proj(z);
    }
    std::cout << w << "\n";
}


void benchStdProj(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::proj(z); 
    }
    std::cout << w << "\n";
}

void benchCpplexPolar(){
    const int N = 1e+9;

    double r = 1.0; 
    double theta = 0.5;
    Complex w = 0;

    for(int i = 0; i < N; i++){
        w += polar(r, theta);
    }
    std::cout << w << "\n";
}

void benchStdPolar(){
    const int N = 1e+9;

    double r = 1.0; 
    double theta = 0.5;
    std::complex<double> w = 0;

    for(int i = 0; i < N; i++){
        w += std::polar(r, theta);
    }
    std::cout << w << "\n";
}

void benchCpplexExp(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += exp(z);
    }
    std::cout << w << "\n";
}

void benchStdExp(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::exp(z); 
    }
    std::cout << w << "\n";
}

void benchCpplexLog(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += log(z);
    }
    std::cout << w << "\n";
}

void benchStdLog(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::log(z);
    }
    std::cout << w << "\n";
}

void benchCpplexLog10(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += log10(z);
    }
    std::cout << w << "\n";
}

void benchStdLog10(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::log10(z);
    }
    std::cout << w << "\n";
}

void benchCpplexPow(){
    const int N = 1e+9;

    Complex z1 = 1 + 1_j;
    Complex z2 = 3 + 3_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += pow(z1, z2);
    }
    std::cout << w << "\n";
}

void benchStdPow(){
    const int N = 1e+9;

    std::complex<double> z1 = 1.0 + 1i;
    std::complex<double> z2 = 3.0 + 3i;

    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::pow(z1, z2); 
    }
    std::cout << w << "\n";
}


void benchCpplexSqrt(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += sqrt(z);
    }
    std::cout << w << "\n";
}

void benchStdSqrt(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;

    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::sqrt(z); 
    }
    std::cout << w << "\n";
}

void benchCpplexSin(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += sin(z);
    }
    std::cout << w << "\n";
}

void benchStdSin(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::sin(z);
    }
    std::cout << w << "\n";
}

void benchCpplexCos(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += cos(z);
    }
    std::cout << w << "\n";
}

void benchStdCos(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::cos(z);
    }
    std::cout << w << "\n";
}

void benchCpplexTan(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += tan(z);
    }
    std::cout << w << "\n";
}

void benchStdTan(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::tan(z);
    }
    std::cout << w << "\n";
}

void benchCpplexAsin(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += asin(z);
    }
    std::cout << w << "\n";
}

void benchStdAsin(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::asin(z); 
    }
    std::cout << w << "\n";
}

void benchCpplexAcos(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += acos(z);
    }
    std::cout << w << "\n";
}

void benchStdAcos(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::acos(z); 
    }
    std::cout << w << "\n";
}

void benchCpplexAtan(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += atan(z);
    }
    std::cout << w << "\n";
}

void benchStdAtan(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::atan(z); 
    }
    std::cout << w << "\n";
}

void benchCpplexSinh(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += sinh(z);
    }
    std::cout << w << "\n";
}

void benchStdSinh(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::sinh(z);
    }
    std::cout << w << "\n";
}

void benchCpplexCosh(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += cosh(z);
    }
    std::cout << w << "\n";
}

void benchStdCosh(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::cosh(z);
    }
    std::cout << w << "\n";
}

void benchCpplexTanh(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += tanh(z);
    }
    std::cout << w << "\n";
}

void benchStdTanh(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::tanh(z);
    }
    std::cout << w << "\n";
}

void benchCpplexAsinh(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += asinh(z);
    }
    std::cout << w << "\n";
}

void benchStdAsinh(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::asinh(z); 
    }
    std::cout << w << "\n";
}

void benchCpplexAcosh(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += acosh(z);
    }
    std::cout << w << "\n";
}

void benchStdAcosh(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::acosh(z); 
    }
    std::cout << w << "\n";
}

void benchCpplexAtanh(){
    const int N = 1e+9;

    Complex z = 1 + 1_j;
    Complex w = 0; 
    for(int i = 0; i < N; i++){
        w += atanh(z);
    }
    std::cout << w << "\n";
}

void benchStdAtanh(){
    const int N = 1e+9;

    std::complex<double> z = 1.0 + 1i;
    std::complex<double> w = 0; 
    for(int i = 0; i < N; i++){
        w += std::atanh(z); 
    }
    std::cout << w << "\n";
}