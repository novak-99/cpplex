#ifndef EXPONENTIAL_HPP
#define EXPONENTIAL_HPP

#include "Complex/Complex.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Exponential : public ContinuousDistribution {
        public:

            constexpr Exponential(const Complex& lambda) noexcept : lambda(lambda) {

            }

            constexpr Complex getLambda() const noexcept {
                return lambda; 
            }

            constexpr void setLambda(const Complex& lambda) noexcept {
                this->lambda = lambda;
            }

            constexpr double pdf(const Complex& z) const noexcept {
                return pdf(z.getX(), lambda.getX()) * pdf(z.getY(), lambda.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), lambda.getX()) * cdf(z.getY(), lambda.getY()); 
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
            
                std::exponential_distribution dX{lambda.getX()};
                std::exponential_distribution dY{lambda.getY()};

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
            inline double pdf(const double x, const double lambda) const noexcept {
                return lambda * std::exp(-lambda * x);
            }

            inline double cdf(const double x, const double lambda) const noexcept {
                return 1 - std::exp(-lambda * x);
            }
            
            inline double entropy(const double lambda) const noexcept {
                return 1 - std::log(lambda);
            }

            Complex lambda; 
    }; 
}


#endif // EXPONENTIAL_HPP