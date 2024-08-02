#ifndef RAYLEIGH_HPP
#define RAYLEIGH_HPP

#include "Complex/Complex.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Rayleigh : public ContinuousDistribution {
        public:

            constexpr Rayleigh(const Complex& sigma) noexcept : sigma(sigma) {

            }

            constexpr Complex getLambda() const noexcept {
                return sigma; 
            }

            constexpr void setLambda(const Complex& sigma) noexcept {
                this->sigma = sigma;
            }

            constexpr double pdf(const Complex& z) const noexcept {
                return pdf(z.getX(), sigma.getX()) * pdf(z.getY(), sigma.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), sigma.getX()) * cdf(z.getY(), sigma.getY()); 
            }

            // inline double logpdf(const Complex& z) const noexcept {
            //     return std::log(pdf(z)); 
            // }

            // inline double logcdf(const Complex& z) const noexcept {
            //     return std::log(cdf(z));  
            // }

            inline std::vector<Complex> rand(const int numSamples) const noexcept {
                std::random_device rd{};
                std::mt19937 gen{rd()};
            
                std::uniform_real_distribution dX{0.0, 1.0};
                std::uniform_real_distribution dY{0.0, 1.0};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = invcdf(dX(gen),  sigma.getX()) + invcdf(dY(gen), sigma.getY()) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(sigma.getX()) + entropy(sigma.getY()); 
            }
        private:
            inline double pdf(const double x, const double sigma) const noexcept {
                return (x/(sigma*sigma)) * std::exp(-x*x / (2 * sigma * sigma)); 
            }

            inline double cdf(const double x, const double sigma) const noexcept {
                return 1 - std::exp(-x*x / (2 * sigma * sigma)); 
            }
            
            inline double entropy(const double sigma) const noexcept {
                return 1 + std::log(sigma / std::sqrt(2)) + eulerGamma.real() / 2; 
            }

            // amazing
            // https://www.math.wm.edu/~leemis/chart/UDR/PDFs/RayleighV.pdf
            inline double invcdf(const double x, const double sigma) const noexcept {
                return std::sqrt(sigma * std::log(1 / (1 - sigma)));
            }

            Complex sigma; 
    }; 
}


#endif // RAYLEIGH_HPP