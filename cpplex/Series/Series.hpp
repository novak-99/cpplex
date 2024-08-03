// #ifndef SERIES_HPP
// #define SERIES_HPP

// #include <cmath>
// #include "Complex.hpp"
// #include "NumericalAnalysis.hpp"
// #include "Constants.hpp"

// namespace cpplex{

//     enum ConvergenceCondition {
//         inconclusive,
//         divergent, 
//         condConvergent, 
//         absConvergent, 
//         inconclusiveLimit 
//     };

//     class Series{

//         public:
//             constexpr Series(Complex (*f)(Complex), Complex lowerLimit, Complex upperLimit) noexcept 
//             : f(f), lowerLimit(lowerLimit), upperLimit(upperLimit) {

//             }

//             constexpr void setF(Complex (*f)(Complex)) noexcept {
//                 this->f = f; 
//             }

//             constexpr void setLowerLimit(const Complex lowerLimit) noexcept {
//                 this->lowerLimit = lowerLimit; 
//             }

//             constexpr void setUpperLimit(const Complex upperLimit) noexcept {
//                 this->upperLimit = upperLimit; 
//             }

//             constexpr Complex (*getF())(Complex) {
//                 return f; 
//             }

//             constexpr Complex getLowerLimit() const noexcept {
//                 return lowerLimit; 
//             }

//             constexpr Complex getUpperLimit() const noexcept {
//                 return upperLimit; 
//             }

//             constexpr ConvergenceCondition divergenceTest() const noexcept {
//                 double eps = 1e-10; 
//                 return mod(limit(f, INF)) < eps ? inconclusive : divergent; 
//             }

//            constexpr ConvergenceCondition ratioTest() const noexcept {
//                 auto ratioFn = [](Complex z, Complex(*f)(Complex) ) { return f(z+1)/f(z); };
//                 double L = mod(limit(ratioFn, f, INF));
//                 std::cout << L << "\n";
//                 if(std::isnan(L)){
//                     return inconclusiveLimit; 
//                 }
//                 else if(L < 1){
//                     return absConvergent; 
//                 }
//                 else if(L > 1){
//                     return divergent; 
//                 }
//                 else {
//                     return inconclusive; 
//                 }
//                 return inconclusive;
//             }

//             constexpr ConvergenceCondition integralTest() const noexcept {
//                 Complex I = integral(f, 1, INF); 
//                 if(isnan(I)){
//                     return divergent; 
//                 }
//                 else{
//                     return absConvergent;
//                 }
//             }

//            constexpr ConvergenceCondition rootTest() const noexcept {
//                 auto ratioFn = [](Complex z, Complex(*f)(Complex) ) { return  pow(mod(f(z)), 1.0/z); };
//                 double L = limit(ratioFn, f, INF).real(); 
//                 std::cout << L << "\n";
//                 if(std::isnan(L)){
//                     return inconclusiveLimit; 
//                 }
//                 else if(L < 1){
//                     return absConvergent; 
//                 }
//                 else if(L > 1){
//                     return divergent; 
//                 }
//                 else {
//                     return inconclusive; 
//                 }
//                 return inconclusive;
//             }



//         private:
//             Complex (*f)(Complex); 
//             Complex lowerLimit; 
//             Complex upperLimit; 

        
//     };
// }

// #endif // SERIES_HPP