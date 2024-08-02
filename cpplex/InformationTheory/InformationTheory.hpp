#ifndef INFORMATIONTHEORY_HPP
#define INFORMATIONTHEORY_HPP

#include "Complex/Complex.hpp"
#include "Constants/Constants.hpp"
#include "FFT/FFT.hpp"
#include "NumericalAnalysis/NumericalAnalysis.hpp"

// ASSUMES P(X, Y) FOR INDEPENDENT X & Y. 
// P(X, Y) = P(X)P(Y). 
// H[X, Y] = H[X] + H[Y] for independent X & Y. 

namespace cpplex{

    inline namespace {

        const double eps = 1e-8; // for log

        constexpr double entropy_(const std::vector<double>& X) noexcept {
            int n = X.size(); 
            double sum = 0; 
            
            for(int i = 0; i < n; i++){
                if(X[i] > 0) {
                    sum += -X[i] * std::log(X[i]); 
                }
                else if(X[i] == 0) ; // lol 
                else return NINF.real(); 
            }
            return sum; 
        }

        constexpr double klDiv_(const std::vector<double>& X, const std::vector<double>& Y) noexcept {
            int n = X.size(); 
            double sum = 0; 

            for(int i = 0; i < n; i++){
                if(X[i] > 0 && Y[i] > 0){
                    sum += X[i] * std::log(X[i] / Y[i]); 
                }
                else if(X[i] == 0 && Y[i] >= 0) ; // lol 
                else return INF.real(); 
            }
            return sum; 
        }

        constexpr double jsDiv_(const std::vector<double>& X, const std::vector<double>& Y) noexcept {
            int n = X.size(); 

            std::vector<double> M(n); 

            for(int i = 0; i < n; i++){
                M[i] = 0.5 * (X[i] + Y[i]); 
            }
            return 0.5 * (klDiv_(X, M) + klDiv_(Y, M)); 
        }
    
        constexpr double shannonInformation_(const double x) noexcept {
            return x > 0 ? -std::log(x) : NINF.real(); 
        }

        constexpr double crossEntropy_(const std::vector<double>& X, const std::vector<double>& Y) noexcept {
            int n = X.size(); 
            double sum = 0; 

            for(int i = 0; i < n; i++){
                if(X[i] < 0 || Y[i] < 0){
                    return INF.real(); 
                }
                else if(X[i] > 0 && Y[i] > 0){
                    sum += -X[i] * std::log(Y[i]); 
                }
            }
            return sum; 
        }
    }

    constexpr std::vector<double> probNorm(const std::vector<double>& X){
        int n = X.size(); 
        std::vector<double> normX(X); 
        double sum = 0; 
        for(int i = 0; i < n; i++){
            if(X[i] < 0) return X; // Do NOTHING. Can't normalize.
            sum += X[i]; 
        }
        if(sum == 0) return X; // Do NOTHING. Can't normalize.
        for(int i = 0; i < n; i++){
            normX[i] /= sum; 
        }
        return normX; 
    }

    constexpr double entropy(const std::vector<Complex>& X) noexcept {
        auto [reX, imX] = decomp(X);
        return entropy_(probNorm(reX)) + entropy_(probNorm(imX)); 
    }

    constexpr double klDiv(const std::vector<Complex>& X, const std::vector<Complex>& Y) noexcept {
        auto [reX, imX] = decomp(X);
        auto [reY, imY] = decomp(Y);
        return klDiv_(probNorm(reX), probNorm(reY)) + klDiv_(probNorm(imX), probNorm(imY)); 
    }

    constexpr double jsDiv(const std::vector<Complex>& X, const std::vector<Complex>& Y) noexcept {
        auto [reX, imX] = decomp(X);
        auto [reY, imY] = decomp(Y);
        return jsDiv_(probNorm(reX), probNorm(reY)) + jsDiv_(probNorm(imX), probNorm(imY)); 
    }

    constexpr double shannonInformation(const Complex& z) noexcept {
        return shannonInformation_(z.real()) + shannonInformation_(z.im()); 
    }

    constexpr double crossEntropy(const std::vector<Complex>& X, const std::vector<Complex>& Y) noexcept {
        auto [reX, imX] = decomp(X);
        auto [reY, imY] = decomp(Y);
        return crossEntropy_(probNorm(reX), probNorm(reY)) + crossEntropy_(probNorm(imX), probNorm(imY)); 
    }

    inline namespace {
        constexpr double entropy_(Complex (*f)(Complex), double a, double b) noexcept {
            auto fe = [](Complex z, std::vector<Complex (*)(Complex)> fargs) { 
                Complex fz = fargs[0](z); 
                return fz * log(fz + eps); // exp(-x^2) can't handle very much leeway 
            }; 
            return integral(fe, a, b, {f}).real();  // MUST BE REAL. 
        }

        constexpr double klDiv_(Complex (*f)(Complex), Complex (*g)(Complex), double a, double b) noexcept {
            auto fe = [](Complex z, std::vector<Complex (*)(Complex)> fargs) { 
                Complex fz = fargs[0](z); 
                return fz * log(fz / (fargs[1](z)+ eps) + eps); 
            }; 
            return integral(fe, a, b, {f, g}).real();  // MUST BE REAL. 
        }

        constexpr double jsDiv_(Complex (*f)(Complex), Complex (*g)(Complex), double a, double b) noexcept {
            auto fe1 = [](Complex z, std::vector<Complex (*)(Complex)> fargs) { 
                Complex fz = fargs[0](z); 
                Complex mz = 0.5 * (fz + fargs[1](z)); 
                return fz * log(fz / (mz + eps) + eps); 
            };

            auto fe2 = [](Complex z, std::vector<Complex (*)(Complex)> fargs) { 
                Complex gz = fargs[1](z); 
                Complex mz = 0.5 * (fargs[0](z) + gz); 
                return gz * log(gz / (mz + eps) + eps); 
            };
            return 0.5 * (integral(fe1, a, b, {f, g}).real() + integral(fe2, a, b, {f, g}).real()); 
        }

        constexpr double crossEntropy_(Complex (*f)(Complex), Complex (*g)(Complex), double a, double b) noexcept {
            auto fe = [](Complex z, std::vector<Complex (*)(Complex)> fargs) { return -fargs[0](z) * log(fargs[1](z) + eps); }; // fe, for f_entropy. 
            return integral(fe, a, b, {f, g}).real();  // MUST BE REAL. 
        }

    }

    constexpr double entropy(Complex (*fr)(Complex), Complex (*fi)(Complex), double a, double b) noexcept {
        return entropy_(fr, a, b) + entropy_(fi, a, b); // probably they have the same support.
    }

    constexpr double klDiv(Complex (*fr)(Complex), Complex (*fi)(Complex), Complex (*gr)(Complex), Complex (*gi)(Complex), double a, double b) noexcept {
        return klDiv_(fr, gr, a, b) + klDiv_(fi, gi, a, b); // probably they have the same support.
    }

    constexpr double jsDiv(Complex (*fr)(Complex), Complex (*fi)(Complex), Complex (*gr)(Complex), Complex (*gi)(Complex), double a, double b) noexcept {
        return jsDiv_(fr, gr, a, b) + jsDiv_(fi, gi, a, b); // probably they have the same support.
    }

    constexpr double crossEntropy(Complex (*fr)(Complex), Complex (*fi)(Complex), Complex (*gr)(Complex), Complex (*gi)(Complex), double a, double b) noexcept {
        return jsDiv_(fr, gr, a, b) + jsDiv_(fi, gi, a, b); // probably they have the same support.
    }

    inline double binaryEntropy(const double p) noexcept {
        return -p * std::log(p + eps) - (1 - p) * std::log(1 - p + eps); 
    }
}

#endif // INFORMATIONTHEORY_HPP