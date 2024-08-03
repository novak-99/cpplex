#ifndef LOGISTIC_HPP
#define LOGISTIC_HPP

#include "Complex.hpp"
#include "Special.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Logistic : public ContinuousDistribution {
        public:

            constexpr Logistic(const Complex& mu, const Complex& s) noexcept : mu(mu), s(s) {

            }

            constexpr Complex getMu() const noexcept {
                return mu; 
            }

            constexpr Complex getS() const noexcept {
                return s; 
            }

            constexpr void setMu(const Complex& mu) noexcept {
                this->mu = mu;
            }

            constexpr void setS(const Complex& s) noexcept {
                this->s = s; 
            }

            // This does assume independence, which isn't correct. I will have to make an updated version which uses the matrix forms for cov., etc. 
            constexpr double pdf(const Complex& z) const noexcept {
                return pdf(z.getX(), mu.getX(), s.getX()) * pdf(z.getY(), mu.getY(), s.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), mu.getX(), s.getX()) * cdf(z.getY(), mu.getY(), s.getY()); 
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
            
                std::uniform_real_distribution dX{0.0, 1.0};
                std::uniform_real_distribution dY{0.0, 1.0};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = invcdf(dX(gen),  mu.getX(), s.getX()) + invcdf(dY(gen), mu.getY(), s.getY()) * 1_j;
                }
                return samples; 
            }
            inline double entropy() const noexcept {
                return entropy(s.getX()) + entropy(s.getY()); 
            }
        private:
            inline double pdf(const double x, const double mu, const double s) const noexcept {
                double sigmoidX = sigmoid((x - mu)/s).real(); 
                return (std::exp(-(x - mu) / s) / s) * sigmoidX * sigmoidX;
            }

            inline double cdf(const double x, const double mu, const double s) const noexcept {
                return sigmoid((x - mu)/s).real(); 
            }
            
            inline double entropy(const double s) const noexcept {
                return 0.5 * std::log(2 * M_PI * s * s) + 0.5; 
            }

            inline double invcdf(const double p, const double mu, const double s) const noexcept {
                return mu + s * std::log(p / (1 - p)); 
            }

            Complex mu; 
            Complex s; 
    }; 
}


#endif // LOGISTIC_HPP