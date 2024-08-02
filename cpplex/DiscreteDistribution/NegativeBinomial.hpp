#ifndef NEGATIVEBINOMIAL_HPP
#define NEGATIVEBINOMIAL_HPP

#include "Complex/Complex.hpp"
#include "Special/Special.hpp"
#include "InformationTheory/InformationTheory.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class NegativeBinomial : public DiscreteDistribution {
        public:

            constexpr NegativeBinomial(const Complex& r, const Complex& p) noexcept : r(r), p(p) {

            }

            constexpr Complex getR() const noexcept {
                return r; 
            }

            constexpr void setR(const Complex& r) noexcept {
                this->r = r; 
            }

            constexpr Complex getP() const noexcept {
                return p; 
            }

            constexpr void setP(const Complex& p) noexcept {
                this->p = p;
            }

            constexpr double pmf(const Complex& z) const noexcept {
                return pmf(z.getX(), r.getX(), p.getX()) * pmf(z.getY(), r.getY(), p.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), r.getY(), p.getX()) * cdf(z.getY(), r.getY(), p.getY()); 
            }

            inline std::vector<Complex> rand(const int numSamples) const noexcept {
                std::random_device rd{};
                std::mt19937 gen{rd()};
            
                std::negative_binomial_distribution dX{int(r.getX()), p.getX()};
                std::negative_binomial_distribution dY{int(r.getY()), p.getY()};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = dX(gen) + dY(gen) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(r.getX(), p.getX()) + entropy(r.getY(), p.getY()); 
            }
        private:
            inline double pmf(const double k, const double r, const double p) const noexcept {
                return binom(k + r - 1, k) * std::pow(1 - p, k) * std::pow((p), r); 
            }

            inline double cdf(const double k, const double r, const double p) const noexcept {
                return betaincReg(p, r, k + 1); 
            }
            
            // https://arxiv.org/pdf/1708.06394
            inline double entropy(const double r, const double p) const noexcept {
                auto f = [](Complex t, std::vector<Complex> args) { 
                    return (pow(1 - t, args[0] - 1) - 1) * (pow((1 + args[1] * t) / (1 - args[1]), args[0]) \
                        + args[0] * args[1] * t / (1 - args[1]) - 1) / (t * log(1 - t)); 
                };
                return r * (binaryEntropy(p) - p * std::log(r)) / (1 - p) + integral(f, 0, 1, {r, p}).real();
            }

            Complex r; 
            Complex p; 
    }; 
}


#endif // NEGATIVEBINOMIAL_HPP