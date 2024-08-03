#include <iostream>
#include <Complex/Complex.hpp>
#include <complex>
#include "NumericalAnalysis/NumericalAnalysis.hpp"
#include "FFT/FFT.hpp"
#include "InformationTheory/InformationTheory.hpp"
#include "Special/Special.hpp"
#include "ContinuousDistribution/ContinuousDistribution.hpp"
#include "ContinuousDistribution/Normal.hpp"
#include "ContinuousDistribution/Gamma.hpp"
#include "ContinuousDistribution/Laplace.hpp"
#include "ContinuousDistribution/Cauchy.hpp"
#include "ContinuousDistribution/LogNormal.hpp"
#include "ContinuousDistribution/Weibull.hpp"
#include "ContinuousDistribution/ChiSquared.hpp"
#include "ContinuousDistribution/Exponential.hpp"
#include "ContinuousDistribution/Uniform.hpp"
#include "ContinuousDistribution/Rayleigh.hpp"
#include "ContinuousDistribution/Logistic.hpp"
#include "ContinuousDistribution/Pareto.hpp"
#include "ContinuousDistribution/Chi.hpp"
#include "ContinuousDistribution/Triangular.hpp"
#include "Transforms/Transforms.hpp"
#include "Signal/Signal.hpp"
#include "DiscreteDistribution/DiscreteDistribution.hpp"
#include "DiscreteDistribution/Poisson.hpp"
#include "DiscreteDistribution/Bernoulli.hpp"
#include "DiscreteDistribution/Binomial.hpp"
#include "DiscreteDistribution/Geometric.hpp"
#include "DiscreteDistribution/NegativeBinomial.hpp"
#include <Benchmarks/Complex.hpp>
#include <Tests/Complex.hpp>
#include <Benchmarks/Zeta.hpp>

using namespace cpplex; 
using namespace std::complex_literals;

// Complex f(Complex z){
//     return exp(-z*z);
// }

// I use this file for running formal unit tests, benchmarks, and for quick and dirty testing 

int main(){
    // benchCpplexZetaSum();

    // std::cout << gamma(1 + 1i) << " " << gamma(1i) << "\n";
    // benchStdAtanh();


        auto fRe = [](Complex t) { return exp(-t * t); }; // Gaussian function.
        auto fIm = [](Complex t) { return exp(-t * t); }; // Gaussian function.
        auto gRe = [](Complex t) { return exp(-t * t); }; // Gaussian function. 
        auto gIm = [](Complex t) { return exp(-t * t); }; // Gaussian function. 
        std::cout << klDiv(fRe, fIm, gRe, gIm, NINF.real(), INF.real()) << "\n"; // Should be ~ 0 (epsilon value included in logs may influence precision).

    return 0; 
    // std::cout << arg(Complex(0.4, -0.0)) << "\n";

    // std::cout << arg(Complex(-0.5, -0.0)) << "\n";

    // std::cout << std::arg(-0.5 + 0.0i) << "\n";


    // std::vector<Complex> a = {1,2,3,4i};
    // std::vector<Complex> b = {1,2,3,4i,5};
    // auto c = crosscorr(a, b); 
    // for(int i = 0; i < c.size(); i++){
    //     std::cout << c[i] << "\n"; 
    // }


    // std::cout << hyp2f1(1, 2, 3, 0.5) << "\n";

    // std::vector<Complex> a{0.5, 0.5}; 
    // std::cout << entropy(a) << "\n";
    // std::cout << yv(0.5, 1i + 1) << "\n";
    // std::cout << conv(fn, fn, 1) << "\n";
    // std::cout << secantMethod(fn, 0.01 + 1i, 0.1 + 1i, 1000) << "\n";

    // std::cout << alpha(1, 1) << "\n";

    // std::cout << weightomega(5 + 5i) << "\n";

    // std::cout << dawsn(1 + 1i) << "\n";

    // std::vector<Complex> X = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}; 
    // std::vector<Complex> X = {2,4,6,8,10};
    // std::vector<Complex> Y = fft(ifft(X));
    // for(int i = 0; i < Y.size(); i++){
    //     std::cout << Y[i] << "\n";
    // }


    // std::vector<Complex> X = {1,2,3,4 + 5_j}; 
    // std::vector<Complex> Y = entropy(X);
    // for(int i = 0; i < Y.size(); i++){
    //     std::cout <<  Y[i] << "\n";
    // }

    

    // std::cout << (1 + 3_j).std() << "\n";
    // std::cout << eulersMethod(cexp, 2, 1, 2_j, 0.000001) << "\n"; 

    // auto roots = rootsOfUnity(1, 5); 

    // for(auto root : roots){
    //     std::cout << root << "\n";
    // }


    // std::complex<double> z_cpp = 3.5 + 5i; 
    // std::complex<double> z_cpp_empty(0,0);

    // Complex z_marc_melikyan = 3.5 + 5_j; 


    // Complex z_empty(0, 0); 
    // for(int i = 0; i < 2 * N; i++){
    //     z_empty = z_marc_melikyan + z_marc_melikyan; 

    // }

    // std::cout << z_empty << "\n";


    return 0; 
}