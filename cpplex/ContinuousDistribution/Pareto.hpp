#ifndef PARETO_HPP
#define PARETO_HPP

#include "Complex/Complex.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Pareto : public ContinuousDistribution {
        public:

            constexpr Pareto(const Complex& xm, const Complex& alpha) noexcept : xm(xm), alpha(alpha) {

            }

            constexpr Complex getXm() const noexcept {
                return xm; 
            }

            constexpr Complex getAlpha() const noexcept {
                return alpha; 
            }

            constexpr void setXm(const Complex& xm) noexcept {
                this->xm = xm;
            }

            constexpr void setAlpha(const Complex& alpha) noexcept {
                this->alpha = alpha; 
            }

            // This does assume independence, which isn't correct. I will have to make an updated version which uses the matrix forms for cov., etc. 
            constexpr double pdf(const Complex& z) const noexcept {
                return pdf(z.getX(), xm.getX(), alpha.getX()) * pdf(z.getY(), xm.getY(), alpha.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), xm.getX(), alpha.getX()) * cdf(z.getY(), xm.getY(), alpha.getY()); 
            }

            inline std::vector<Complex> rand(const int numSamples) const noexcept {
                std::random_device rd{};
                std::mt19937 gen{rd()};
            
                std::uniform_real_distribution dX{0.0, 1.0};
                std::uniform_real_distribution dY{0.0, 1.0};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = invcdf(dX(gen),  xm.getX(), alpha.getX()) + invcdf(dY(gen), xm.getY(), alpha.getY()) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(xm.getX(), alpha.getX()) + entropy(xm.getY(), alpha.getY()); 
            }
        private:
            inline double pdf(const double x, const double xm, const double alpha) const noexcept {
                return alpha * std::pow(xm, alpha) / std::pow(x, alpha + 1); 
            }

            inline double cdf(const double x, const double xm, const double alpha) const noexcept {
                return 1 - std::pow(xm / x, alpha); 
            }
            
            inline double entropy(const double xm, const double alpha) const noexcept {
                return std::log(xm / alpha + std::exp(1 + 1 / alpha)); 
            }

            inline double invcdf(const double p, const double xm, const double alpha) const noexcept {
                return alpha * std::pow(p, alpha - 1) * xm / alpha; 
            }

            Complex xm; 
            Complex alpha; 
    }; 
}


#endif // PARETO_HPP