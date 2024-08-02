#ifndef TRIANGULAR_HPP
#define TRIANGULAR_HPP

#include "Complex/Complex.hpp"
#include <vector>
#include <cmath>
#include <random>

namespace cpplex{
    class Triangular : public ContinuousDistribution {
        public:

            constexpr Triangular(const Complex& a, const Complex& b, const Complex& c) noexcept : a(a), b(b), c(c) {

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
                return pdf(z.getX(), a.getX(), b.getX(), c.getX()) * pdf(z.getY(), a.getY(), b.getY(), c.getY()); 
            }

            inline double cdf(const Complex& z) const noexcept {
                return cdf(z.getX(), a.getX(), b.getX(), c.getY()) * cdf(z.getY(), a.getY(), b.getY(), c.getY()); 
            }

            inline std::vector<Complex> rand(const int numSamples) const noexcept {
                std::random_device rd{};
                std::mt19937 gen{rd()};
            
                std::uniform_real_distribution dX{0.0, 1.0};
                std::uniform_real_distribution dY{0.0, 1.0};

                std::vector<Complex> samples(numSamples); 
                for(int i = 0; i < numSamples; i++){
                    samples[i] = invcdf(dX(gen), a.getX(), b.getX(), c.getX()) + invcdf(dY(gen), a.getY(), b.getY(), c.getY()) * 1_j;
                }
                return samples; 
            }

            inline double entropy() const noexcept {
                return entropy(a.getX(), b.getX()) + entropy(a.getY(), b.getY()); 
            }
        private:
            inline double pdf(const double x, const double a, const double b, const double c) const noexcept {
                if(x < a) return 0; 
                else if(x >= a && x < c){
                    return 2 * (x - a) / ((b - a) * (c - a)); 
                }
                else if(x == c) return 2 / (b - a); 
                else if(x > c && x <= b){
                    return 2 * (b - x) / ((b - a) * (b - c)); 
                }
                else return 0; 
            }

            inline double cdf(const double x, const double a, const double b, const double c) const noexcept {
                if(x <= a) return 0; 
                else if(x > a && x <= c){
                    return 2 * (x - a) * (x - a) / ((b - a) * (b - c)); 
                }
                else if(x > c && x < b){
                    return 1 - ((b - x) * (b - x)) / ((b - a) * (b - c)); 
                }
                else return 1;  
            }
            
            inline double entropy(const double a, const double b) const noexcept {
                return 0.5 + std::log((b - a) / 2); 
            }

            // https://www.drdawnwright.com/easy-excel-inverse-triangular-distribution-for-monte-carlo-simulations/
            inline double invcdf(const double p, const double a, const double b, const double c) const noexcept {
                if(p <= (c - a) / (b - a)) return a + std::sqrt(p * (b - a) * (c - a));
                else return b - std::sqrt((1 -p) * (b - a) * (b - c));
            }

            Complex a; 
            Complex b; 
            Complex c; 
    }; 
}


#endif // TRIANGULAR_HPP