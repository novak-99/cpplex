#ifndef LAPLACE_HPP
#define LAPLACE_HPP

#include "Complex.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Laplace : public ContinuousDistribution {
        public:

            constexpr Laplace(const Complex& mu, const Complex& b) noexcept : mu(mu), b(b) {

            }

            constexpr Complex getMu() const noexcept {
                return mu; 
            }

            constexpr Complex getB() const noexcept {
                return b; 
            }

            constexpr void setMu(const Complex& mu) noexcept {
                this->mu = mu;
            }

            constexpr void setB(const Complex& b) noexcept {
                this->b = b; 
            }

            // This does assume independence, which isn't correct. I will have to make an updated version which uses the matrix forms for cov., etc. 
            constexpr double pdf(const Complex& z) const noexcept {
                return pdf(z.getX(), mu.getX(), b.getX()) * pdf(z.getY(), mu.getY(), b.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), mu.getX(), b.getX()) * cdf(z.getY(), mu.getY(), b.getY()); 
            }

            inline std::vector<Complex> rand(const int numSamples) const noexcept {
                std::random_device rd{};
                std::mt19937 gen{rd()};
            
                std::uniform_real_distribution dX{0.0, 1.0};
                std::uniform_real_distribution dY{0.0, 1.0};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = invcdf(dX(gen),  mu.getX(), b.getX()) + invcdf(dY(gen), mu.getY(), b.getY()) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(b.getX()) + entropy(b.getY()); 
            }
        private:
            inline double pdf(const double x, const double mu, const double b) const noexcept {
                double expTerm(-std::abs(x - mu) / b); 
                return expTerm / (2 * b); 
            }

            inline double cdf(const double x, const double mu, const double b) const noexcept {
                if(x < mu){
                    return 0.5 * std::exp((x - mu) / b); 
                } 
                else if(x == mu){
                    return 0.5; 
                }
                else{
                    return 1 - 0.5 * std::exp(-(x - mu) / b); 
                }
            }
            
            inline double entropy(const double b) const noexcept {
                return std::log(2 * b) + 1; 
            }

            inline double invcdf(const double p, const double mu, const double b) const noexcept {
                return mu - b * sgn(p - 0.5) * std::log(1 - 2 * std::abs(p - 0.5)); 
            }

            Complex mu; 
            Complex b; 
    }; 
}


#endif // LAPLACE_HPP