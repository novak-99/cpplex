#ifndef NORMAL_HPP
#define NORMAL_HPP

#include "Complex.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Normal : public ContinuousDistribution {
        public:

            constexpr Normal(const Complex& mu, const Complex& sigma) noexcept : mu(mu), sigma(sigma) {

            }

            constexpr Complex getMu() const noexcept {
                return mu; 
            }

            constexpr Complex getSigma() const noexcept {
                return sigma; 
            }

            constexpr void setMu(const Complex& mu) noexcept {
                this->mu = mu;
            }

            constexpr void setSigma(const Complex& sigma) noexcept {
                this->sigma = sigma; 
            }

            // This does assume independence, which isn't correct. I will have to make an updated version which uses the matrix forms for cov., etc. 
            constexpr double pdf(const Complex& z) const noexcept {
                return pdf(z.getX(), mu.getX(), sigma.getX()) * pdf(z.getY(), mu.getY(), sigma.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), mu.getX(), sigma.getX()) * cdf(z.getY(), mu.getY(), sigma.getY()); 
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
            
                std::normal_distribution dX{mu.getX(), sigma.getX()};
                std::normal_distribution dY{mu.getY(), sigma.getY()};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = dX(gen) + dY(gen) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(sigma.getX()) + entropy(sigma.getY()); 
            }
        private:
            inline double pdf(const double x, const double mu, const double sigma) const noexcept {
                double expTerm(-0.5 * std::pow((std::log(x) - mu)/sigma, 2)); 
                return expTerm / (sigma * std::sqrt(2 * M_PI) * x); 
            }

            inline double cdf(const double x, const double mu, const double sigma) const noexcept {
                return 0.5 * ((1 + std::erf(std::log(x) - mu) / (std::sqrt(2) * sigma))); 
            }
            
            inline double entropy(const double mu, const double sigma) const noexcept {
                return mu + 0.5 * std::log(2 * M_PI * std::exp(1.0) * sigma * sigma); 
            }

            Complex mu; 
            Complex sigma; 
    }; 
}


#endif // NORMAL_HPP