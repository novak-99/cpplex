#ifndef CHI_HPP
#define CHI_HPP

#include "Complex/Complex.hpp"
#include "Special/Special.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Chi : public ContinuousDistribution {
        public:

            constexpr Chi(const Complex& k) noexcept : k(k) {

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
                    samples[i] = std::sqrt(dX(gen)) + std::sqrt(dY(gen)) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(k.getX()) + entropy(k.getY()); 
            }
        private:
            inline double pdf(const double x, const double k) const noexcept {
                return (std::pow(x, k - 1) * std::exp(-x * x/2)) / (std::pow(2, k/2 - 1) * gamma(k/2).real()); 
            }

            inline double cdf(const double x, const double k) const noexcept {
                return pdf(x * x / 2, k / 2); 
            }
            
            inline double entropy(const double k) const noexcept {
                return std::log(gamma(k / 2).real()) + 0.5 * (k - std::log(2) - (k - 1) * psi(k / 2).real()); 
            }

            Complex k; 
    }; 
}


#endif // CHI_HPP