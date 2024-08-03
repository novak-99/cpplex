#ifndef BERNOULLI_HPP
#define BERNOULLI_HPP

#include "Complex.hpp"
#include "Special.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Bernoulli : public DiscreteDistribution {
        public:

            constexpr Bernoulli(const Complex& p) noexcept : p(p) {

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
            
                std::bernoulli_distribution dX{p.getX()};
                std::bernoulli_distribution dY{p.getY()};

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
                if(k == 0) return 1 - p;  
                else if(k == 1) return p; 
                else return 0; 
            }

            inline double cdf(const double k, const double p) const noexcept {
                if(k < 0) return 0;  
                else if(k > 0 && k <= 1) return 1 - p; 
                else return 1;  
            }
            
            inline double entropy(const double p) const noexcept {
                return -(1 - p) * std::log(1 - p) - p * std::log(p); 
            }

            Complex p; 
    }; 
}


#endif // BERNOULLI_HPP