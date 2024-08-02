#ifndef UNIFORM_HPP
#define UNIFORM_HPP

#include "Complex/Complex.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex {
    class Uniform : public ContinuousDistribution {
        public:

            constexpr Uniform(const Complex& a, const Complex& b) noexcept : a(a), b(b) {

            }

            constexpr Complex getA() const noexcept {
                return a; 
            }

            constexpr Complex getB() const noexcept {
                return b; 
            }

            constexpr void setA(const Complex& a) noexcept {
                this->a = a;
            }

            constexpr void setB(const Complex& b) noexcept {
                this->b = b; 
            }

            // This does assume independence, which isn't correct. I will have to make an updated version which uses the matrix forms for cov., etc. 
            constexpr double pdf(const Complex& z) const noexcept {
                return pdf(z.getX(), a.getX(), b.getX()) * pdf(z.getY(), a.getY(), b.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), a.getX(), b.getX()) * cdf(z.getY(), a.getY(), b.getY()); 
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
            
                std::uniform_real_distribution dX{a.getX(), b.getX()};
                std::uniform_real_distribution dY{a.getY(), b.getY()};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = dX(gen) + dY(gen) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(a.getX(), b.getX()) + entropy(a.getX(), b.getY()); 
            }
        private:
            inline double pdf(const double x, const double a, const double b) const noexcept {
                return (x >= a && x <= b) ?  1 / (b-a) : 0; 
            }

            inline double cdf(const double x, const double a, const double b) const noexcept {
                if(x < a) return 0; 
                else if(x >= a && x <= b) return (x - a) / (b - a); 
                else return 1; 
            }
            
            inline double entropy(const double a, const double b) const noexcept {
                return std::log(b - a); 
            }

            Complex a; 
            Complex b; 
    }; 
}


#endif // UNIFORM_HPP