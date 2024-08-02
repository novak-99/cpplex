#ifndef GEOMETRIC_HPP
#define GEOMETRIC_HPP

#include "Complex/Complex.hpp"
#include "Special/Special.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Geometric : public DiscreteDistribution {
        public:

            constexpr Geometric(const Complex& p) noexcept : p(p) {

            }

            constexpr Complex getLambda() const noexcept {
                return p; 
            }

            constexpr void setLambda(const Complex& p) noexcept {
                this->p = p;
            }

            constexpr double pmf(const Complex& z) const noexcept {
                return pmf(z.getX(), p.getX()) * pmf(z.getY(), p.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), p.getX()) * cdf(z.getY(), p.getY()); 
            }

            inline std::vector<Complex> rand(const int numSamples) const noexcept {
                std::random_device rd{};
                std::mt19937 gen{rd()};
            
                std::geometric_distribution dX{p.getX()};
                std::geometric_distribution dY{p.getY()};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = dX(gen) + dY(gen) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(p.getX()) + entropy(p.getY()); 
            }
        private:
            inline double pmf(const double k, const double p) const noexcept {
                return std::pow((1 - p), k - 1) * p; 
            }

            inline double cdf(const double k, const double p) const noexcept {
                return k >= 1 ? 1 - std::pow(1 - p, std::floor(k)) : 0; 
            }
            
            inline double entropy(const double p) const noexcept {
                return (-(1 - p) * std::log(1 - p) - p * std::log(p)) / p; 
            }

            Complex p; 
    }; 
}


#endif // GEOMETRIC_HPP