#ifndef WEIBULL_HPP
#define WEIBULL_HPP

#include "Complex.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Weibull : public ContinuousDistribution {
        public:

            constexpr Weibull(const Complex& lambda, const Complex& k) noexcept : lambda(lambda), k(k) {

            }

            constexpr Complex getLambda() const noexcept {
                return lambda; 
            }

            constexpr Complex getK() const noexcept {
                return k; 
            }

            constexpr void setLambda(const Complex& lambda) noexcept {
                this->lambda = lambda;
            }

            constexpr void setK(const Complex& k) noexcept {
                this->k = k; 
            }

            constexpr double pdf(const Complex& z) const noexcept {
                return pdf(z.getX(), lambda.getX(), k.getX()) * pdf(z.getY(), lambda.getY(), k.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), lambda.getX(), k.getX()) * cdf(z.getY(), lambda.getY(), k.getY()); 
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
            
                std::weibull_distribution dX{lambda.getX(), k.getX()};
                std::weibull_distribution dY{lambda.getY(), k.getY()};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = dX(gen) + dY(gen) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(lambda.getX(), k.getX()) + entropy(lambda.getY(), k.getY()); 
            }
        private:
            inline double pdf(const double x, const double lambda, const double k) const noexcept {
                return x < 0 ? 0 : (k/lambda) * std::pow(x/lambda, k - 1) * std::exp(std::pow(-x / lambda, k)); 
            }

            inline double cdf(const double x, const double lambda, const double k) const noexcept {
                return x < 0 ? 0 : 1 - std::exp(std::pow(-x / lambda, k)); 
            }
            
            inline double entropy(const double lambda, const double k) const noexcept {
                return eulerGamma.real() * (1 - 1/k) + std::log(lambda / k) + 1; 
            }

            Complex lambda; 
            Complex k; 
    }; 
}


#endif // WEIBULL_HPP