#ifndef CONTINUOUSDISTRIBUTION_HPP
#define CONTINUOUSDISTRIBUTION_HPP

#include "Complex/Complex.hpp"

#include <vector>

// Please note: this distribution class assumes the Re and Im components are INDEPENDENT. 
namespace cpplex{
    class ContinuousDistribution {
        public:
            // REQUIRED. 
            virtual constexpr double pdf(const Complex& z) const noexcept = 0; 
            virtual constexpr double cdf(const Complex& z) const noexcept = 0; 

            inline double logpdf(const Complex& z) const noexcept {
                return std::log(pdf(z)); 
            }

            inline double logcdf(const Complex& z) const noexcept {
                return std::log(cdf(z));  
            }
            virtual constexpr std::vector<Complex> rand(const int numSamples) const noexcept = 0; 
            virtual constexpr double entropy() const noexcept = 0; 
        private: 
    }; 
}


#endif // CONTINUOUSDISTRIBUTION_HPP