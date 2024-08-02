#include <Complex/Complex.hpp>
#include <Constants/Constants.hpp>
#include <iostream>
#include <cassert>
#include <cmath>

using namespace cpplex;


void testCompoundAdditionComplex(){
    std::vector<std::pair<Complex, Complex>> z = { {2, 2}, {2, 2 + 4_j}, {2 + 4_j, 2 + 4_j} }; 
    std::vector<Complex> compoundAdditionZ = { 4, 4 + 4_j, 4 + 8_j };

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        Complex zf = z[i].first; 
        zf += z[i].second; 
        assert(zf == compoundAdditionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testCompoundAdditionReal(){
    std::vector<std::pair<Complex, double>> z = { {2, 2}, {2 + 4_j, 2}, {1_j, 2} }; 
    std::vector<Complex> compoundAdditionZ = { 4, 4 + 4_j, 2 + 1_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        Complex zf = z[i].first; 
        zf += z[i].second; 
        assert(zf == compoundAdditionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testCompoundSubtractionComplex(){
    std::vector<std::pair<Complex, Complex>> z = { {2, 2}, {2, 2 + 4_j}, {2 + 4_j, 2 + 4_j} }; 
    std::vector<Complex> compoundSubtractionZ = { 0, -4_j, 0 };

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        Complex zf = z[i].first; 
        zf -= z[i].second; 
        assert(zf == compoundSubtractionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testCompoundSubtractionReal(){
    std::vector<std::pair<Complex, double>> z = { {2, 2}, {2 + 4_j, 2}, {1_j, 2} }; 
    std::vector<Complex> compoundSubtractionZ = { 0, 4_j, -2 + 1_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        Complex zf = z[i].first; 
        zf -= z[i].second; 
        assert(zf == compoundSubtractionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testCompoundMultiplicationComplex(){
    std::vector<std::pair<Complex, Complex>> z = { {2, 2}, {2, 2 + 4_j}, {2 + 4_j, 2 + 4_j} }; 
    std::vector<Complex> compoundMultiplicationZ = { 4, 4 + 8_j, -12 + 16_j };

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        Complex zf = z[i].first; 
        zf *= z[i].second; 
        assert(zf == compoundMultiplicationZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testCompoundMultiplicationReal(){
    std::vector<std::pair<Complex, double>> z = { {2, 2}, {2 + 4_j, 2}, {1_j, 2} }; 
    std::vector<Complex> compoundMultiplicationZ = { 4, 4 + 8_j, 2_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        Complex zf = z[i].first; 
        zf *= z[i].second; 
        assert(zf == compoundMultiplicationZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testCompoundDivisionComplex(){
    std::vector<std::pair<Complex, Complex>> z = { {2, 2}, {2, 2 + 4_j}, {2 + 4_j, 3 + 1_j} }; 
    std::vector<Complex> compoundDivisionZ = { 1, 0.2 - 0.4_j, 1 + 1_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        Complex zf = z[i].first; 
        zf /= z[i].second; 
        assert(zf == compoundDivisionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testCompoundDivisionReal(){
    std::vector<std::pair<Complex, double>> z = { {2, 2}, {2 + 4_j, 2}, {1_j, 2} }; 
    std::vector<Complex> compoundDivisionZ = { 1, 1 + 2_j, 0.5_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        Complex zf = z[i].first; 
        zf /= z[i].second; 
        assert(zf == compoundDivisionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testAdditionComplexComplex(){
    std::vector<std::pair<Complex, Complex>> z = { {2, 2}, {2, 2 + 4_j}, {2 + 4_j, 2 + 4_j} }; 
    std::vector<Complex> additionZ = { 4, 4 + 4_j, 4 + 8_j };

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(z[i].first + z[i].second == additionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testAdditionComplexReal(){
    std::vector<std::pair<Complex, double>> z = { {2, 2}, {2 + 4_j, 2}, {4_j, 2} }; 
    std::vector<Complex> additionZ = { 4, 4 + 4_j, 2 + 4_j };

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(z[i].first + z[i].second == additionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testAdditionRealComplex(){
    std::vector<std::pair<double, Complex>> z = { {2, 2}, {2, 2 + 4_j}, {2, 4_j} }; 
    std::vector<Complex> additionZ = { 4, 4 + 4_j, 2 + 4_j };

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(z[i].first + z[i].second == additionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testSubtractionComplexComplex(){
    std::vector<std::pair<Complex, Complex>> z = { {2, 2}, {2, 2 + 4_j}, {2 + 4_j, 2 + 4_j} }; 
    std::vector<Complex> subtractionZ = { 0, -4_j, 0 };

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(z[i].first - z[i].second == subtractionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testSubtractionComplexReal(){
    std::vector<std::pair<Complex, double>> z = { {2, 2}, {2 + 4_j, 2}, {4_j, 2} }; 
    std::vector<Complex> subtractionZ = { 0, 4_j, -2 + 4_j };

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(z[i].first - z[i].second == subtractionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testSubtractionRealComplex(){
    std::vector<std::pair<double, Complex>> z = { {2, 2}, {2, 2 + 4_j}, {2, 4_j} }; 
    std::vector<Complex> subtractionZ = { 0, -4_j, 2 - 4_j };

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(z[i].first - z[i].second == subtractionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testMultiplicationComplexComplex(){
    std::vector<std::pair<Complex, Complex>> z = { {2, 2}, {2, 2 + 4_j}, {2 + 4_j, 2 + 4_j} }; 
    std::vector<Complex> multiplicationZ = { 4, 4 + 8_j, -12 + 16_j };

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(z[i].first * z[i].second == multiplicationZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testMultiplicationComplexReal(){
    std::vector<std::pair<Complex, double>> z = { {2, 2}, {2 + 4_j, 2}, {4_j, 2} }; 
    std::vector<Complex> multiplicationZ = { 4, 4 + 8_j, 8_j };

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(z[i].first * z[i].second == multiplicationZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testMultiplicationRealComplex(){
    std::vector<std::pair<double, Complex>> z = { {2, 2}, {2, 2 + 4_j}, {2, 4_j} }; 
    std::vector<Complex> multiplicationZ = { 4, 4 + 8_j, 8_j };

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(z[i].first * z[i].second == multiplicationZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testDivisionComplexComplex(){
    std::vector<std::pair<Complex, Complex>> z = { {2, 2}, {2, 2 + 4_j}, {2 + 4_j, 3 + 1_j} }; 
    std::vector<Complex> divisionZ = { 1, 0.2 - 0.4_j, 1 + 1_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(z[i].first / z[i].second == divisionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testDivisionComplexReal(){
    std::vector<std::pair<Complex, double>> z = { {2, 2}, {2 + 4_j, 2}, {4_j, 2} }; 
    std::vector<Complex> divisionZ = { 1, 1 + 2_j, 2_j };

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(z[i].first / z[i].second == divisionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testDivisionRealComplex(){
    std::vector<std::pair<double, Complex>> z = { {2, 2}, {2, 2 + 4_j}, {2, 4_j} }; 
    std::vector<Complex> divisionZ = { 1, 0.2 - 0.4_j, -0.5_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(z[i].first / z[i].second == divisionZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testMod(){
    std::vector<Complex> z = {3 + 5_j, 0.5 + 0.55_j, 1 + 50_j, -3 - 4_j, -0.5 + 0.001_j}; 
    std::vector<double> modZ = {std::sqrt(34), std::sqrt(0.5525), std::sqrt(2501), 5, std::sqrt(0.250001)};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(mod(z[i]) == modZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testArg(){
    std::vector<Complex> z = {1_j, 0, 2 + 2_j, -1 -1_j, -1_j}; 
    std::vector<double> argZ = {M_PI/2, 0, M_PI/4, -3 * M_PI/4, -M_PI/2};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(arg(z[i]) == argZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testNorm(){
    std::vector<Complex> z = {3 + 5_j, 0.5 + 0.55_j, 1 + 50_j, -3 - 4_j, -0.5 + 0.001_j}; 
    std::vector<double> normZ = {34, 0.5525, 2501, 25, 0.250001};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(norm(z[i]) == normZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testConj(){
    std::vector<Complex> z = {3 + 5_j, 0.5 + 0.55_j, 1 + 50_j, -3 - 4_j, -0.5 + 0.001_j}; 
    std::vector<Complex> conjZ = {3 - 5_j, 0.5 - 0.55_j, 1 - 50_j, -3 + 4_j, -0.5 - 0.001_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(conj(z[i]) == conjZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testDot(){
    std::vector<std::pair<Complex, Complex>> z = { {3 + 5_j, 0}, {1_j, 1_j}, {1 + 1_j, 2}, {1 + 1_j, 1 + 1_j}, {0, 0}}; 
    std::vector<double> dotZ = {0, 1, 2, 2, 0};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(dot(z[i].first, z[i].second) == dotZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testProj(){
    std::vector<Complex> z = {0, 4 + 4_j, INF, NINF, Complex(0, NINF.real())}; 
    std::vector<Complex> projZ = {0, 4 + 4_j, INF, INF, Complex(INF.real(), -0.0)};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(proj(z[i]) == projZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testPolar(){
    const double ERR = 1e-8;
    std::vector<std::pair<double, double>> z = {  {1, 0}, {2, M_PI}, {3, M_PI/2} }; 
    std::vector<Complex> polarZ = {1, -2, 3_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = polar(z[i].first, z[i].second).real(); 
        double y = polar(z[i].first, z[i].second).im(); 

        assert(std::abs(x - polarZ[i].getX()) < ERR && std::abs(y - polarZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

// maybe use a different test case.
void testExp(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {1, M_PI * 1_j}; 
    std::vector<Complex> expZ = {std::exp(1), -1};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = exp(z[i]).real();
        double y = exp(z[i]).im();

        assert(std::abs(x - expZ[i].getX()) < ERR && std::abs(y - expZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testLog(){
    std::vector<Complex> z = {std::exp(1), 1_j}; 
    std::vector<Complex> logZ = {1, M_PI/2 * 1_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(log(z[i]) == logZ[i]);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testLog10(){
    std::vector<Complex> z = {std::exp(1), 1_j}; 
    std::vector<Complex> logZ = {1, M_PI/2 * 1_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(log10(z[i]) == logZ[i] / std::log(10));
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testLog2(){
    std::vector<Complex> z = {std::exp(1), 1_j}; 
    std::vector<Complex> logZ = {1, M_PI/2 * 1_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        assert(log2(z[i]) == logZ[i] / std::log(2));
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testSqrt(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {0, 25, -1, 2_j, 3 + 4_j }; 
    std::vector<Complex> sqrtZ = {0, 5, 1_j, 1 + 1_j, 2 + 1_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = sqrt(z[i]).real();
        double y = sqrt(z[i]).im();

        assert(std::abs(x - sqrtZ[i].getX()) < ERR && std::abs(y - sqrtZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testPowComplexComplex(){
    const double ERR = 1e-8;
    std::vector<std::pair<Complex, Complex>> z = { {5, 2}, {std::exp(1), log(1_j)}, {1_j, 1_j}}; 
    std::vector<Complex> powZ = {25, 1_j, std::exp(-M_PI/2)};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = pow(z[i].first, z[i].second).real();
        double y = pow(z[i].first, z[i].second).im();

        // assert(std::abs(x - sqrtZ[i].getX()) < ERR && std::abs(y - sqrtZ[i].getY()) < ERR);

        assert(std::abs(x - powZ[i].getX()) < ERR && std::abs(y - powZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testPowRealComplex(){
    const double ERR = 1e-8;
    std::vector<std::pair<double, Complex>> z = { {5,2}, {std::exp(1.0), 1_j}, {std::exp(1), log(1_j)} }; // pow(real, complex)
    std::vector<Complex> powZ = {25, std::cos(1) + std::sin(1) * 1_j, 1_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = pow(z[i].first, z[i].second).real();
        double y = pow(z[i].first, z[i].second).im();

        // assert(std::abs(x - sqrtZ[i].getX()) < ERR && std::abs(y - sqrtZ[i].getY()) < ERR);

        assert(std::abs(x - powZ[i].getX()) < ERR && std::abs(y - powZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testPowComplexReal(){
    const double ERR = 1e-8;
    std::vector<std::pair<Complex, double>> z = { {5,2}, {1_j, 2}, {1 + 1_j, 2} }; // pow(complex, real)
    std::vector<Complex> powZ = {25, -1, 2_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = pow(z[i].first, z[i].second).real();
        double y = pow(z[i].first, z[i].second).im();

        // assert(std::abs(x - sqrtZ[i].getX()) < ERR && std::abs(y - sqrtZ[i].getY()) < ERR);

        assert(std::abs(x - powZ[i].getX()) < ERR && std::abs(y - powZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testSin(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {0, 1_j * std::log(1 + std::sqrt(2)),  1_j * M_PI,  1_j * std::log(2 + std::sqrt(3)),  M_PI/2 + 1_j * std::log(1 + std::sqrt(2))}; 
    std::vector<Complex> sinZ = {0, 1_j, 1_j * std::sinh(M_PI), std::sqrt(3) * 1_j, std::sqrt(2)};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = sin(z[i]).real();
        double y = sin(z[i]).im();

        assert(std::abs(x - sinZ[i].getX()) < ERR && std::abs(y - sinZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testCos(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {0, 1_j * std::log(1 + std::sqrt(2)),  1_j * M_PI, 1_j * std::log(2 + std::sqrt(3)),  M_PI/2 + 1_j * std::log(1 + std::sqrt(2))}; 
    std::vector<Complex> cosZ = {1, std::sqrt(2), std::cosh(M_PI), 2, -1_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = cos(z[i]).real();
        double y = cos(z[i]).im();

        assert(std::abs(x - cosZ[i].getX()) < ERR && std::abs(y - cosZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testTan(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {0, 1_j * std::log(1 + std::sqrt(2)),  1_j * M_PI,  1_j * std::log(2 + std::sqrt(3)),  M_PI/2 + 1_j * std::log(1 + std::sqrt(2))}; 
    std::vector<Complex> tanZ = {0, std::sqrt(2)/2 * 1_j, std::tanh(M_PI) * 1_j, std::sqrt(3)/2 * 1_j, std::sqrt(2) * 1_j};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = tan(z[i]).real();
        double y = tan(z[i]).im();

        assert(std::abs(x - tanZ[i].getX()) < ERR && std::abs(y - tanZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testAsin(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {0, 1_j, 1_j * std::sinh(M_PI), std::sqrt(3) * 1_j, std::sqrt(2)};
    std::vector<Complex> asinZ = {0, 1_j * std::log(1 + std::sqrt(2)),  1_j * M_PI,  1_j * std::log(2 + std::sqrt(3)),  M_PI/2 + 1_j * std::log(1 + std::sqrt(2))}; 

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = asin(z[i]).real();
        double y = asin(z[i]).im();

        assert(std::abs(x - asinZ[i].getX()) < ERR && std::abs(y - asinZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testAcos(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {1, std::sqrt(2), std::cosh(M_PI), 2, -1_j};
    std::vector<Complex> acosZ = {0, -1_j * std::log(1 + std::sqrt(2)),  -1_j * M_PI, -1_j * std::log(2 + std::sqrt(3)),  M_PI/2 + 1_j * std::log(1 + std::sqrt(2))}; 

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = acos(z[i]).real();
        double y = acos(z[i]).im();

        assert(std::abs(x - acosZ[i].getX()) < ERR && std::abs(y - acosZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testAtan(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {0, std::sqrt(2)/2 * 1_j, std::tanh(M_PI) * 1_j, std::sqrt(3)/2 * 1_j, std::sqrt(2) * 1_j};
    std::vector<Complex> atanZ = {0, 1_j * std::log(1 + std::sqrt(2)),  1_j * M_PI,  1_j * std::log(2 + std::sqrt(3)),  M_PI/2 + 1_j * std::log(1 + std::sqrt(2))}; 

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = atan(z[i]).real();
        double y = atan(z[i]).im();

        assert(std::abs(x - atanZ[i].getX()) < ERR && std::abs(y - atanZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testSinh(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {0, log(1_j),  1_j * M_PI,  std::log(2) + 1_j * M_PI/2,  std::log(2 + std::sqrt(3))}; 
    std::vector<Complex> sinhZ = {0, 1_j, 0, 1.25_j, std::sqrt(3)};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = sinh(z[i]).real();
        double y = sinh(z[i]).im();

        assert(std::abs(x - sinhZ[i].getX()) < ERR && std::abs(y - sinhZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testCosh(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {0, log(1_j),  1_j * M_PI,  std::log(2) + 1_j * M_PI/2,  std::log(2 + std::sqrt(3))}; 
    std::vector<Complex> coshZ = {1, 0, -1, 0.75_j, 2};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = cosh(z[i]).real();
        double y = cosh(z[i]).im();

        assert(std::abs(x - coshZ[i].getX()) < ERR && std::abs(y - coshZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testTanh(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {0, log(2_j),  1_j * M_PI,  std::log(2) + 1_j * M_PI/2,  std::log(2 + std::sqrt(3))}; 
    std::vector<Complex> tanhZ = {0, 1.66666667, 0, 1.66666667, std::sqrt(3)/2};

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = tanh(z[i]).real();
        double y = tanh(z[i]).im();

        assert(std::abs(x - tanhZ[i].getX()) < ERR && std::abs(y - tanhZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testAsinh(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {0, 1_j, std::sqrt(3), 1.25_j, std::sqrt(3)};
    std::vector<Complex> asinhZ = {0, log(1_j),  std::log(2 + std::sqrt(3)),  std::log(2) + 1_j * M_PI/2,  std::log(2 + std::sqrt(3))}; 

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = asinh(z[i]).real();
        double y = asinh(z[i]).im();

        assert(std::abs(x - asinhZ[i].getX()) < ERR && std::abs(y - asinhZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testAcosh(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {1, 0, -1, 0.75_j, 2};
    std::vector<Complex> acoshZ = {0, log(1_j),  1_j * M_PI,  std::log(2) + 1_j * M_PI/2,  std::log(2 + std::sqrt(3))}; 

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = acosh(z[i]).real();
        double y = acosh(z[i]).im();

        assert(std::abs(x - acoshZ[i].getX()) < ERR && std::abs(y - acoshZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}

void testAtanh(){
    const double ERR = 1e-8;
    std::vector<Complex> z = {0, 1.66666667, 1_j, 1.66666667, std::sqrt(3)/2};
    std::vector<Complex> atanhZ = {0, log(2_j), 1_j * M_PI/4,  std::log(2) + 1_j * M_PI/2,  std::log(2 + std::sqrt(3))}; 

    int n = z.size(); 
    for(int i = 0; i < n; i++){
        double x = atanh(z[i]).real();
        double y = atanh(z[i]).im();

        assert(std::abs(x - atanhZ[i].getX()) < ERR && std::abs(y - atanhZ[i].getY()) < ERR);
        std::cout << "Passed Test " << i + 1 << "\n";
    }
}