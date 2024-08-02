#ifndef MAXWELLBOLTZMAN_HPP
#define MAXWELLBOLTZMAN_HPP

#include "Complex/Complex.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class MaxwellBoltzman : public ContinuousDistribution {
        public:

            constexpr MaxwellBoltzman(const Complex& a) noexcept : a(a) {

            }

            constexpr Complex getA() const noexcept {
                return a; 
            }

            constexpr void setA(const Complex& a) noexcept {
                this->a = a;
            }

            constexpr double pdf(const Complex& z) const noexcept {
                return pdf(z.getX(), a.getX()) * pdf(z.getY(), a.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), a.getX()) * cdf(z.getY(), a.getY()); 
            }

            // inline double logpdf(const Complex& z) const noexcept {
            //     return std::log(pdf(z)); 
            // }

            // inline double logcdf(const Complex& z) const noexcept {
            //     return std::log(cdf(z));  
            // }

            inline std::vector<Complex> rand(const int numSamples) const noexcept {

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = rand(a.getX()) + rand(a.getY()) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(a.getX()) + entropy(a.getY()); 
            }
        private:
            inline double pdf(const double x, const double a) const noexcept {
                double x2 = x * x; 
                double a2 = a * a; 
                return std::sqrt(2 / M_PI) * (x2 / (a2 * a)) * std::exp(-x2 / (2 * a2)); 
            }

            inline double cdf(const double x, const double a) const noexcept {
                return std::erf(x / (std::sqrt(2) * a)) - std::sqrt(2 / M_PI) * (x / a) * std::exp(-x * x / (2 * a * a)); 
            }
            
            inline double entropy(const double a) const noexcept {
                return std::log(a * std::sqrt(2 * M_PI)) + eulerGamma.real() - 0.5; 
            }

            inline double rand(const double a) const noexcept {
                std::random_device rd{};
                std::mt19937 gen{rd()};

                std::maxwellboltzman_distribution d{0, a.getX()};

                double vx = d(gen);
                double vy = d(gen);
                double vz = d(gen);

                return std::sqrt(vx * vx + vy * vy + vz * vz); 
            }

            Complex a; 
    }; 
}


#endif // MAXWELLBOLTZMAN_HPP