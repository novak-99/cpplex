#ifndef CHISQUARED_HPP
#define CHISQUARED_HPP

#include "Complex.hpp"
#include "Special.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class ChiSquared : public ContinuousDistribution {
        public:

            constexpr ChiSquared(const Complex& k) noexcept : k(k) {

            }

            constexpr Complex getK() const noexcept {
                return k; 
            }

            constexpr void setK(const Complex& k) noexcept {
                this->k = k;
            }

            constexpr double pdf(const Complex& z) const noexcept {
                return pdf(z.getX(), k.getX()) * pdf(z.getY(), k.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), k.getX()) * cdf(z.getY(), k.getY()); 
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
            
                std::chi_squared_distribution dX{k.getX()};
                std::chi_squared_distribution dY{k.getY()};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = dX(gen) + dY(gen) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(k.getX()) + entropy(k.getY()); 
            }
        private:
            inline double pdf(const double x, const double k) const noexcept {
                return (std::pow(2, k/2) * std::pow(x, k/2 - 1) * std::exp(-x/2)) / gamma(k/2).real(); 
            }

            inline double cdf(const double x, const double k) const noexcept {
                return gammainc(k/2, x/2) / gamma(k/2).real(); 
            }
            
            inline double entropy(const double k) const noexcept {
                return k/2 - std::log(2 * gamma(k/2).real()) + (1 - k/2) * psi(k/2).real(); 
            }

            Complex k; 
    }; 
}


#endif // CHISQUARED_HPP