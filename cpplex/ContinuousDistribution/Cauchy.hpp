#ifndef CAUCHY_HPP
#define CAUCHY_HPP

#include "Complex/Complex.hpp"
#include "Special/Special.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Cauchy : public ContinuousDistribution {
        public:

            constexpr Cauchy(const Complex& mu, const Complex& gamma) noexcept : mu(mu), gamma(gamma) {

            }

            constexpr Complex getMu() const noexcept {
                return mu; 
            }

            constexpr Complex setGamma() const noexcept {
                return gamma; 
            }

            constexpr void setMu(const Complex& mu) noexcept {
                this->mu = mu;
            }

            constexpr void setGamma(const Complex& gamma) noexcept {
                this->gamma = gamma; 
            }

            constexpr double pdf(const Complex& z) const noexcept {
                return pdf(z.getX(), mu.getX(), gamma.getX()) * pdf(z.getY(), mu.getY(), gamma.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), mu.getX(), gamma.getX()) * cdf(z.getY(), mu.getY(), gamma.getY()); 
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
            
                std::cauchy_distribution dX{mu.getX(), gamma.getX()};
                std::cauchy_distribution dY{mu.getY(), gamma.getY()};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = dX(gen) + dY(gen) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(mu.getX(), gamma.getX()) + entropy(mu.getY(), gamma.getY()); 
            }
        private:
            inline double pdf(const double x, const double mu, const double gamma) const noexcept {
                double normX = (x - mu) / gamma; 
                return 1 / ( (M_PI * gamma) * (1 + normX * normX) ); 
            }

            inline double cdf(const double x, const double mu, const double gamma) const noexcept {
                double normX = (x - mu) / gamma; 
                return (1 / M_PI) * std::atan(normX) + 0.5; 
            }
            
            inline double entropy(const double mu, const double gamma) const noexcept {
                return std::log(4 * M_PI * gamma); 
            }

            Complex mu; 
            Complex gamma; 
    }; 
}


#endif // CAUCHY_HPP