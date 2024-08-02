#ifndef GAMMA_HPP
#define GAMMA_HPP

#include "Complex/Complex.hpp"
#include "Special/Special.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Gamma : public ContinuousDistribution {
        public:

            constexpr Gamma(const Complex& alpha, const Complex& beta) noexcept : alpha(alpha), beta(beta) {

            }

            constexpr Complex getAlpha() const noexcept {
                return alpha; 
            }

            constexpr Complex getBeta() const noexcept {
                return beta; 
            }

            constexpr void setAlpha(const Complex& alpha) noexcept {
                this->alpha = alpha;
            }

            constexpr void setBeta(const Complex& beta) noexcept {
                this->beta = beta; 
            }

            constexpr double pdf(const Complex& z) const noexcept {
                return pdf(z.getX(), alpha.getX(), beta.getX()) * pdf(z.getY(), alpha.getY(), beta.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), alpha.getX(), beta.getX()) * cdf(z.getY(), alpha.getY(), beta.getY()); 
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
            
                std::gamma_distribution dX{alpha.getX(), beta.getX()};
                std::gamma_distribution dY{alpha.getY(), beta.getY()};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = dX(gen) + dY(gen) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(alpha.getX(), beta.getX()) + entropy(alpha.getY(), beta.getY()); 
            }
        private:
            inline double pdf(const double x, const double alpha, const double beta) const noexcept {
                return (std::pow(beta, alpha) * std::pow(x, alpha - 1) * std::exp(-beta * x)) / gamma(alpha).real(); 
            }

            inline double cdf(const double x, const double alpha, const double beta) const noexcept {
                return gammainc(alpha, beta * x) / gamma(alpha).real(); 
            }
            
            inline double entropy(const double alpha, const double beta) const noexcept {
                return alpha - std::log(beta) + loggamma(alpha).real() + (1 - alpha) * psi(alpha).real(); 
            }

            Complex alpha; 
            Complex beta; 
    }; 
}


#endif // GAMMA_HPP