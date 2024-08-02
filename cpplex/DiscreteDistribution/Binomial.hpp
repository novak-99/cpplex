#ifndef BINOMIAL_HPP
#define BINOMIAL_HPP

#include "Complex/Complex.hpp"
#include "Special/Special.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Binomial : public DiscreteDistribution {
        public:

            constexpr Binomial(const Complex& n, const Complex& p) noexcept : n(n), p(p) {

            }

            constexpr Complex getN() const noexcept {
                return n; 
            }

            constexpr void setN(const Complex& n) noexcept {
                this->n = n; 
            }

            constexpr Complex getP() const noexcept {
                return p; 
            }

            constexpr void setP(const Complex& p) noexcept {
                this->p = p;
            }

            constexpr double pmf(const Complex& z) const noexcept {
                return pmf(z.getX(), n.getX(), p.getX()) * pmf(z.getY(), n.getY(), p.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), n.getY(), p.getX()) * cdf(z.getY(), n.getY(), p.getY()); 
            }

            inline std::vector<Complex> rand(const int numSamples) const noexcept {
                std::random_device rd{};
                std::mt19937 gen{rd()};
            
                std::binomial_distribution dX{int(n.getX()), p.getX()};
                std::binomial_distribution dY{int(n.getY()), p.getY()};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = dX(gen) + dY(gen) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(n.getX(), p.getX()) + entropy(n.getY(), p.getY()); 
            }
        private:
            inline double pmf(const double k, const double n, const double p) const noexcept {
                if(std::round(k) != k) return 0;
                if(k > n) return 0; 
                return binom(n, k) * std::pow(p, k) * std::pow((1 - p), n - k); 
            }

            inline double cdf(const double k, const double n, const double p) const noexcept {
                return betaincReg(1 - p, n - std::floor(k), 1 + std::floor(p));
            }
            
            inline double entropy(const double n, const double p) const noexcept {
                return 0.5 * std::log(2 * M_PI * std::exp(1.0) * n * p * (1 - p)); 
            }

            Complex n; 
            Complex p; 
    }; 
}


#endif // BINOMIAL_HPP