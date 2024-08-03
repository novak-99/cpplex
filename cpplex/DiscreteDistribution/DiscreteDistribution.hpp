#ifndef DISCRETEDISTRIBUTION_HPP
#define DISCRETEDISTRIBUTION_HPP

#include "Complex.hpp"

#include <vector>

// Please note: this distribution class assumes the Re and Im components are INDEPENDENT. 
namespace cpplex{
    class DiscreteDistribution {
        public:
            // REQUIRED. 
            virtual constexpr double pmf(const Complex& z) const noexcept = 0; 
            virtual constexpr double cdf(const Complex& z) const noexcept = 0; 

            inline double logpmf(const Complex& z) const noexcept {
                return std::log(pmf(z)); 
            }

            inline double logcdf(const Complex& z) const noexcept {
                return std::log(cdf(z));  
            }
            virtual constexpr std::vector<Complex> rand(const int numSamples) const noexcept = 0; 
            virtual constexpr double entropy() const noexcept = 0; 
        private: 
    }; 
}


#endif // DISCRETEDISTRIBUTION_HPP