#ifndef POISSON_HPP
#define POISSON_HPP

#include "Complex/Complex.hpp"
#include "Special/Special.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Poisson : public DiscreteDistribution {
        public:

            constexpr Poisson(const Complex& lambda) noexcept : lambda(lambda) {

            }

            constexpr Complex getLambda() const noexcept {
                return lambda; 
            }

            constexpr void setLambda(const Complex& lambda) noexcept {
                this->lambda = lambda;
            }

            constexpr double pmf(const Complex& z) const noexcept {
                return pmf(z.getX(), lambda.getX()) * pmf(z.getY(), lambda.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), lambda.getX()) * cdf(z.getY(), lambda.getY()); 
            }

            inline std::vector<Complex> rand(const int numSamples) const noexcept {
                std::random_device rd{};
                std::mt19937 gen{rd()};
            
                std::poisson_distribution dX{lambda.getX()};
                std::poisson_distribution dY{lambda.getY()};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = dX(gen) + dY(gen) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(lambda.getX()) + entropy(lambda.getY()); 
            }
        private:
            inline double pmf(const double k, const double lambda) const noexcept {
                if(std::round(k) != k) return 0; // INVALID PMF input! 
                return std::pow(lambda, k) * std::exp(-lambda) / factorial(k); 
            }

            inline double cdf(const double k, const double lambda) const noexcept {
                return gammaincc(std::floor(k + 1), lambda) / factorial(std::floor(k)); 
            }
            
            inline double entropy(const double lambda) const noexcept {
                // PLEASE NOTE -- APPROX entropy for Poisson. Err noted. 
                return 0.5 * std::log(2 * M_PI * std::exp(1.0) * lambda) - 1.0 / (12 * lambda) - 1 / (24 * lambda * lambda) \
                - 19 / (lambda * lambda * lambda); // + O(1 / l^4)
            }

            Complex lambda; 
    }; 
}


#endif // POISSON_HPP