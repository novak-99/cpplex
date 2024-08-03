#ifndef NUMERICALANALYSIS_HPP
#define NUMERICALANALYSIS_HPP

#include "Complex.hpp"
#include "Constants.hpp"

namespace cpplex{

    enum IntegralBound {
        proper, 
        twoSidedInf, 
        pinf,
        ninf,
    };

    class NumericalAnalysis{

    }; 

    constexpr Complex derivative(Complex (*f)(Complex), const Complex& z) noexcept {
        const double h = 1e-8; 
        return (f(z + h) - f(z - h)) / (2 * h); 
    }

//    constexpr Complex realDerivative(Complex (*f)(Complex), const Complex& z) noexcept {
//         return derivative(f, z); 
//     }

//     constexpr Complex imDerivative(Complex (*f)(Complex), const Complex& z) noexcept {
//         return 1_j * derivative(f, z); 
//     }

    constexpr Complex secondDerivative(Complex (*f)(Complex), const Complex& z) noexcept {
        const double h = 1e-5; 
        return ( f(z + h) - 2 * f(z) + f(z - h) ) / (h * h); 
    }

    constexpr Complex constantApproximation(Complex(*f)(Complex), Complex z0) noexcept {
        return f(z0);
    }

    constexpr Complex linearApproximation(Complex(*f)(Complex), Complex z0, Complex z) noexcept {
        return constantApproximation(f, z0) + derivative(f, z0) * (z - z0);
    }

    constexpr Complex quadraticApproximation(Complex(*f)(Complex), Complex z0, Complex z) noexcept {
        return linearApproximation(f, z0, z) + 0.5 * secondDerivative(f, z0) * (z - z0) * (z - z0);
    }

    // // These are the unmixed partials.
    // constexpr Complex realSecondDerivative(Complex (*f)(Complex), const Complex& z) noexcept {
    //     return secondDerivative(f, z); 
    // }

    // constexpr Complex imSecondDerivative(Complex (*f)(Complex), const Complex& z) noexcept {
    //     return  -secondDerivative(f, z); 
    // }

    // constexpr Complex laplacian(Complex (*f)(Complex), const Complex& z) noexcept {
    //     // Will always return 0
    //     return realSecondDerivative(f, z) + imSecondDerivative(f, z); 
    // } 

    constexpr Complex gradientDescent(Complex (*f)(Complex), const Complex& z0, const double alpha, const int epochs) noexcept {
        Complex z(z0); 
        for(int i = 0; i < epochs; i++){
            z -= alpha * derivative(f, z); 
        }
        return z; 
    }

    constexpr Complex secantMethod(Complex (*f)(Complex), const Complex& z0, const Complex& z1, const int epochs) noexcept {
        Complex z(z1); 
        Complex zPrev(z0); 
        for(int i = 0; i < epochs; i++){
            Complex fz = f(z); 
            Complex fzPrev = f(zPrev); 
            z -= f(z) * (z - zPrev) / (fz - fzPrev); 
        }
        return z; 
    }

    constexpr Complex newtonsMethod(Complex (*f)(Complex), const Complex& z0, const int epochs) noexcept {
        Complex z(z0); 
        for(int i = 0; i < epochs; i++){
            z -= f(z) / derivative(f, z); 
        }
        return z; 
    }

    constexpr Complex halleysMethod(Complex (*f)(Complex), const Complex& z0, const int epochs) noexcept {
        Complex z(z0); 
        for(int i = 0; i < epochs; i++){
            Complex df = derivative(f, z); 
            Complex d2f = secondDerivative(f, z);
            Complex fz = f(z); 
            Complex den = 2 * df - fz * d2f; 
            z -= (2 * fz * df) / den; 
        }
        return z; 
    }

    // constexpr Complex eulersMethod(Complex (*derivative)(Complex), const Complex& zf, const Complex& z0, const Complex& f0, const double h) noexcept {

    //     Complex f = f0;  
    //     Complex z = zf;    
    //     while(z != zf){
    //         f += h * derivative(z) * (1 + 1_j); 
    //         z += h;

    //     }

    //     return f; 
    // }

    // rewrite these 3-4 with helper functions. 

    constexpr Complex integral(Complex (*f)(Complex), const Complex& a, const Complex& b) noexcept {
        
        std::vector<double> kronrodNodes = {
            0.991455371120813, 
            -0.991455371120813,
            0.949107912342759,
            -0.949107912342759,
            0.864864423359769,
            -0.864864423359769,
            0.741531185599394,
            -0.741531185599394,
            0.586087235467691,
            -0.586087235467691,
            0.405845151377397,
            -0.405845151377397,
            0.207784955007898,
            -0.207784955007898,
            0.000000000000000
        }; 

        std::vector<double> weights = {
            0.022935322010529, 
            0.022935322010529,
            0.063092092629979,
            0.063092092629979,
            0.104790010322250,
            0.104790010322250,
            0.140653259715525,
            0.140653259715525,
            0.169004726639267,
            0.169004726639267,
            0.190350578064785,
            0.190350578064785,
            0.204432940075298,
            0.204432940075298,
            0.209482141084728,
            0.209482141084728
        }; 

        IntegralBound bounds = proper; 

        Complex af(a); 
        Complex bf(b); 

        if(a == NINF || b == INF){
            af = 0; 
            bf = 1; 
            if(a == NINF && b == INF){
                bounds = twoSidedInf; 
            } 
            else if(b == INF){
                bounds = pinf; 
            }
            else if(a == NINF){
                bounds = ninf; 
            } 
        }

        Complex integral = 0; 
        for(int i = 0; i < kronrodNodes.size(); i++){
            Complex x = ( (kronrodNodes[i] + 1) / 2 ) * (bf - af) + af; 
            Complex w = weights[i] * ( (bf - af) / 2 );

            Complex fx; 
            if(bounds != proper) { // improper integral!
                Complex tau = (1 - x) / x; // substitution for change of integral bounds.. 
                Complex dtau = 1/(x*x); 
                if(bounds == twoSidedInf){
                    fx = (f(tau) + f(-tau)) * dtau; 
                }
                else if(bounds == pinf){
                    fx = f(a + tau) * dtau; 
                }
                else if(bounds == ninf){
                     fx = f(b - tau) * dtau; 
                }
            }
            else{ // else proper integral
                fx = f(x); 
            }

            integral += fx * w; 
        }
        return integral;
    }

    // An overload for functions w/ func
    constexpr Complex integral(Complex (*f)(Complex, std::vector<Complex (*)(Complex)>), const Complex& a, const Complex& b, const std::vector<Complex (*)(Complex)>& fargs) noexcept {
        
        std::vector<double> kronrodNodes = {
            0.991455371120813, 
            -0.991455371120813,
            0.949107912342759,
            -0.949107912342759,
            0.864864423359769,
            -0.864864423359769,
            0.741531185599394,
            -0.741531185599394,
            0.586087235467691,
            -0.586087235467691,
            0.405845151377397,
            -0.405845151377397,
            0.207784955007898,
            -0.207784955007898,
            0.000000000000000
        }; 

        std::vector<double> weights = {
            0.022935322010529, 
            0.022935322010529,
            0.063092092629979,
            0.063092092629979,
            0.104790010322250,
            0.104790010322250,
            0.140653259715525,
            0.140653259715525,
            0.169004726639267,
            0.169004726639267,
            0.190350578064785,
            0.190350578064785,
            0.204432940075298,
            0.204432940075298,
            0.209482141084728,
            0.209482141084728
        }; 

        IntegralBound bounds = proper; 

        Complex af(a); 
        Complex bf(b); 

        if(a == NINF || b == INF){
            af = 0; 
            bf = 1; 
            if(a == NINF && b == INF){
                bounds = twoSidedInf; 
            } 
            else if(b == INF){
                bounds = pinf; 
            }
            else if(a == NINF){
                bounds = ninf; 
            } 
        }

        Complex integral = 0; 
        for(int i = 0; i < kronrodNodes.size(); i++){
            Complex x = ( (kronrodNodes[i] + 1) / 2 ) * (bf - af) + af; 
            Complex w = weights[i] * ( (bf - af) / 2 );

            Complex fx; 
            if(bounds != proper) { // improper integral!
                Complex tau = (1 - x) / x; // substitution for change of integral bounds.. 
                Complex dtau = 1/(x*x); 
                if(bounds == twoSidedInf){
                    fx = (f(tau, fargs) + f(-tau, fargs)) * dtau; 
                }
                else if(bounds == pinf){
                    fx = f(a + tau, fargs) * dtau; 
                }
                else if(bounds == ninf){
                     fx = f(b - tau, fargs) * dtau; 
                }
            }
            else{ // else proper integral
                fx = f(x, fargs); 
            }

            integral += fx * w; 
        }
        return integral;
    }


    // An overload for functions w/ extra args. 
    constexpr Complex integral(Complex (*f)(Complex, std::vector<Complex>), const Complex& a, const Complex& b, const std::vector<Complex>& args) noexcept {
        
        std::vector<double> kronrodNodes = {
            0.991455371120813, 
            -0.991455371120813,
            0.949107912342759,
            -0.949107912342759,
            0.864864423359769,
            -0.864864423359769,
            0.741531185599394,
            -0.741531185599394,
            0.586087235467691,
            -0.586087235467691,
            0.405845151377397,
            -0.405845151377397,
            0.207784955007898,
            -0.207784955007898,
            0.000000000000000
        }; 

        std::vector<double> weights = {
            0.022935322010529, 
            0.022935322010529,
            0.063092092629979,
            0.063092092629979,
            0.104790010322250,
            0.104790010322250,
            0.140653259715525,
            0.140653259715525,
            0.169004726639267,
            0.169004726639267,
            0.190350578064785,
            0.190350578064785,
            0.204432940075298,
            0.204432940075298,
            0.209482141084728,
            0.209482141084728
        }; 

        IntegralBound bounds = proper; 

        Complex af(a); 
        Complex bf(b); 

        if(a == NINF || b == INF){
            af = 0; 
            bf = 1; 
            if(a == NINF && b == INF){
                bounds = twoSidedInf; 
            } 
            else if(b == INF){
                bounds = pinf; 
            }
            else if(a == NINF){
                bounds = ninf; 
            } 
        }

        Complex integral = 0; 
        for(int i = 0; i < kronrodNodes.size(); i++){
            Complex x = ( (kronrodNodes[i] + 1) / 2 ) * (bf - af) + af; 
            Complex w = weights[i] * ( (bf - af) / 2 );

            Complex fx; 
            if(bounds != proper) { // improper integral!
                Complex tau = (1 - x) / x; // substitution for change of integral bounds.. 
                Complex dtau = 1/(x*x); 
                if(bounds == twoSidedInf){
                    fx = (f(tau, args) + f(-tau, args)) * dtau; 
                }
                else if(bounds == pinf){
                    fx = f(a + tau, args) * dtau; 
                }
                else if(bounds == ninf){
                     fx = f(b - tau, args) * dtau; 
                }
            }
            else{ // else proper integral
                fx = f(x, args); 
            }

            integral += fx * w; 
        }
        return integral;
    }

    // An overload for functions w/ extra args + func. 
    constexpr Complex integral(Complex (*f)(Complex,  const std::vector<Complex (*)(Complex)>& fargs, std::vector<Complex>), const Complex& a, const Complex& b,  const std::vector<Complex (*)(Complex)>& fargs, const std::vector<Complex>& args) noexcept {
        
        std::vector<double> kronrodNodes = {
            0.991455371120813, 
            -0.991455371120813,
            0.949107912342759,
            -0.949107912342759,
            0.864864423359769,
            -0.864864423359769,
            0.741531185599394,
            -0.741531185599394,
            0.586087235467691,
            -0.586087235467691,
            0.405845151377397,
            -0.405845151377397,
            0.207784955007898,
            -0.207784955007898,
            0.000000000000000
        }; 

        std::vector<double> weights = {
            0.022935322010529, 
            0.022935322010529,
            0.063092092629979,
            0.063092092629979,
            0.104790010322250,
            0.104790010322250,
            0.140653259715525,
            0.140653259715525,
            0.169004726639267,
            0.169004726639267,
            0.190350578064785,
            0.190350578064785,
            0.204432940075298,
            0.204432940075298,
            0.209482141084728,
            0.209482141084728
        }; 

        IntegralBound bounds = proper; 

        Complex af(a); 
        Complex bf(b); 

        if(a == NINF || b == INF){
            af = 0; 
            bf = 1; 
            if(a == NINF && b == INF){
                bounds = twoSidedInf; 
            } 
            else if(b == INF){
                bounds = pinf; 
            }
            else if(a == NINF){
                bounds = ninf; 
            } 
        }

        Complex integral = 0; 
        for(int i = 0; i < kronrodNodes.size(); i++){
            Complex x = ( (kronrodNodes[i] + 1) / 2 ) * (bf - af) + af; 
            Complex w = weights[i] * ( (bf - af) / 2 );

            Complex fx; 
            if(bounds != proper) { // improper integral!
                Complex tau = (1 - x) / x; // substitution for change of integral bounds.. 
                Complex dtau = 1/(x*x); 
                if(bounds == twoSidedInf){
                    fx = (f(tau, fargs, args) + f(-tau, fargs, args)) * dtau; 
                }
                else if(bounds == pinf){
                    fx = f(a + tau, fargs, args) * dtau; 
                }
                else if(bounds == ninf){
                     fx = f(b - tau, fargs, args) * dtau; 
                }
            }
            else{ // else proper integral
                fx = f(x, fargs, args); 
            }

            integral += fx * w; 
        }
        return integral;
    }

    constexpr Complex integral(Complex (*f)(Complex, std::vector<Complex (*)(Complex)>, std::vector<Complex>), const Complex& a, const Complex& b, const std::vector<Complex (*)(Complex)>& fargs, const std::vector<Complex>& args) noexcept {
        
        std::vector<double> kronrodNodes = {
            0.991455371120813, 
            -0.991455371120813,
            0.949107912342759,
            -0.949107912342759,
            0.864864423359769,
            -0.864864423359769,
            0.741531185599394,
            -0.741531185599394,
            0.586087235467691,
            -0.586087235467691,
            0.405845151377397,
            -0.405845151377397,
            0.207784955007898,
            -0.207784955007898,
            0.000000000000000
        }; 

        std::vector<double> weights = {
            0.022935322010529, 
            0.022935322010529,
            0.063092092629979,
            0.063092092629979,
            0.104790010322250,
            0.104790010322250,
            0.140653259715525,
            0.140653259715525,
            0.169004726639267,
            0.169004726639267,
            0.190350578064785,
            0.190350578064785,
            0.204432940075298,
            0.204432940075298,
            0.209482141084728,
            0.209482141084728
        }; 

        IntegralBound bounds = proper; 

        Complex af(a); 
        Complex bf(b); 

        if(a == NINF || b == INF){
            af = 0; 
            bf = 1; 
            if(a == NINF && b == INF){
                bounds = twoSidedInf; 
            } 
            else if(b == INF){
                bounds = pinf; 
            }
            else if(a == NINF){
                bounds = ninf; 
            } 
        }

        Complex integral = 0; 
        for(int i = 0; i < kronrodNodes.size(); i++){
            Complex x = ( (kronrodNodes[i] + 1) / 2 ) * (bf - af) + af; 
            Complex w = weights[i] * ( (bf - af) / 2 );

            Complex fx; 
            if(bounds != proper) { // improper integral!
                Complex tau = (1 - x) / x; // substitution for change of integral bounds.. 
                Complex dtau = 1/(x*x); 
                if(bounds == twoSidedInf){
                    fx = (f(tau, fargs, args) + f(-tau, fargs, args)) * dtau; 
                }
                else if(bounds == pinf){
                    fx = f(a + tau, fargs, args) * dtau; 
                }
                else if(bounds == ninf){
                     fx = f(b - tau, fargs, args) * dtau; 
                }
            }
            else{ // else proper integral
                fx = f(x, fargs, args); 
            }

            integral += fx * w; 
        }
        return integral;
    }
    

    // An overload for functions w/ extra args + func. 
    constexpr Complex integral(Complex (*f)(Complex, std::vector<Complex (*)(Complex)>, std::vector<Complex>), const Complex& a, const Complex& b, double tau,  const std::vector<Complex (*)(Complex)>& fargs, const std::vector<Complex>& args) noexcept {
        
        std::vector<double> kronrodNodes = {
            0.991455371120813, 
            -0.991455371120813,
            0.949107912342759,
            -0.949107912342759,
            0.864864423359769,
            -0.864864423359769,
            0.741531185599394,
            -0.741531185599394,
            0.586087235467691,
            -0.586087235467691,
            0.405845151377397,
            -0.405845151377397,
            0.207784955007898,
            -0.207784955007898,
            0.000000000000000
        }; 

        std::vector<double> kronrodWeights = {
            0.022935322010529, 
            0.022935322010529,
            0.063092092629979,
            0.063092092629979,
            0.104790010322250,
            0.104790010322250,
            0.140653259715525,
            0.140653259715525,
            0.169004726639267,
            0.169004726639267,
            0.190350578064785,
            0.190350578064785,
            0.204432940075298,
            0.204432940075298,
            0.209482141084728,
            0.209482141084728
        }; 

        std::vector<double> gaussNodes = {
            -0.949107912342759, -0.741531185599394, -0.405845151377397,
            0.000000000000000,  0.405845151377397,  0.741531185599394,  0.949107912342759
        };

        std::vector<double> gaussWeights = {
            0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469,
            0.381830050505119, 0.279705391489277, 0.129484966168870
        };

        IntegralBound bounds = proper; 

        Complex af(a); 
        Complex bf(b); 

        if(a == NINF || b == INF){
            af = 0; 
            bf = 1; 
            if(a == NINF && b == INF){
                bounds = twoSidedInf; 
            } 
            else if(b == INF){
                bounds = pinf; 
            }
            else if(a == NINF){
                bounds = ninf; 
            } 
        }

        bounds = twoSidedInf; 
        // bounds = pinf;

        Complex kronrodIntegral = 0; 
        for(int i = 0; i < kronrodNodes.size(); i++){
            Complex x = ( (kronrodNodes[i] + 1) / 2 ) * (bf - af) + af; 
            Complex w = kronrodWeights[i] * ( (bf - af) / 2 );

            Complex fx; 
            if(bounds != proper) { // improper integral!
                Complex tau = (1 - x) / x; // substitution for change of integral bounds.. 
                Complex dtau = 1/(x*x); 
                if(bounds == twoSidedInf){
                    fx = (f(tau, fargs, args) + f(-tau, fargs, args)) * dtau; 
                }
                else if(bounds == pinf){
                    fx = f(a + tau, fargs, args) * dtau; 
                }
                else if(bounds == ninf){
                     fx = f(b - tau, fargs, args) * dtau; 
                }
            }
            else{ // else proper integral
                fx = f(x, fargs, args); 
            }

            kronrodIntegral += fx * w; 
        }

        Complex gaussIntegral = 0; 
        for(int i = 0; i < gaussNodes.size(); i++){
            Complex x = ( (gaussNodes[i] + 1) / 2 ) * (bf - af) + af; 
            Complex w = gaussWeights[i] * ( (bf - af) / 2 );

            Complex fx; 
            if(bounds != proper) { // improper integral!
                Complex tau = (1 - x) / x; // substitution for change of integral bounds.. 
                Complex dtau = 1/(x*x); 
                if(bounds == twoSidedInf){
                    fx = (f(tau, fargs, args) + f(-tau, fargs, args)) * dtau; 
                }
                else if(bounds == pinf){
                    fx = f(a + tau, fargs, args) * dtau; 
                }
                else if(bounds == ninf){
                     fx = f(b - tau, fargs, args) * dtau; 
                }
            }
            else{ // else proper integral
                fx = f(x, fargs, args); 
            }

            gaussIntegral += fx * w; 
        }

        double epsilon = mod(kronrodIntegral - gaussIntegral); 

        if(epsilon > tau){
            Complex m = (af + bf) / 2; 
            kronrodIntegral = integral(f, af, m, tau/2, fargs, args) + integral(f, m, bf, tau/2, fargs, args); 
        }
        return kronrodIntegral;
    }

// An overload for functions w/ extra args + func. 
    constexpr Complex integral(Complex (*f)(Complex, std::vector<Complex>), const Complex& a, const Complex& b, double tau, const std::vector<Complex>& args) noexcept {
        
        std::vector<double> kronrodNodes = {
            0.991455371120813, 
            -0.991455371120813,
            0.949107912342759,
            -0.949107912342759,
            0.864864423359769,
            -0.864864423359769,
            0.741531185599394,
            -0.741531185599394,
            0.586087235467691,
            -0.586087235467691,
            0.405845151377397,
            -0.405845151377397,
            0.207784955007898,
            -0.207784955007898,
            0.000000000000000
        }; 

        std::vector<double> kronrodWeights = {
            0.022935322010529, 
            0.022935322010529,
            0.063092092629979,
            0.063092092629979,
            0.104790010322250,
            0.104790010322250,
            0.140653259715525,
            0.140653259715525,
            0.169004726639267,
            0.169004726639267,
            0.190350578064785,
            0.190350578064785,
            0.204432940075298,
            0.204432940075298,
            0.209482141084728,
            0.209482141084728
        }; 

        std::vector<double> gaussNodes = {
            -0.949107912342759, -0.741531185599394, -0.405845151377397,
            0.000000000000000,  0.405845151377397,  0.741531185599394,  0.949107912342759
        };

        std::vector<double> gaussWeights = {
            0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469,
            0.381830050505119, 0.279705391489277, 0.129484966168870
        };

        IntegralBound bounds = proper; 

        Complex af(a); 
        Complex bf(b); 

        if(a == NINF || b == INF){
            af = 0; 
            bf = 1; 
            if(a == NINF && b == INF){
                bounds = twoSidedInf; 
            } 
            else if(b == INF){
                bounds = pinf; 
            }
            else if(a == NINF){
                bounds = ninf; 
            } 
        }

        // bounds = twoSidedInf; 
        // bounds = pinf;
        // bounds = proper;

        Complex kronrodIntegral = 0; 
        for(int i = 0; i < kronrodNodes.size(); i++){
            Complex x = ( (kronrodNodes[i] + 1) / 2 ) * (bf - af) + af; 
            Complex w = kronrodWeights[i] * ( (bf - af) / 2 );

            Complex fx; 
            if(bounds != proper) { // improper integral!
                Complex tau = (1 - x) / x; // substitution for change of integral bounds.. 
                Complex dtau = 1/(x*x); 
                if(bounds == twoSidedInf){
                    fx = (f(tau, args) + f(-tau, args)) * dtau; 
                }
                else if(bounds == pinf){
                    fx = f(a + tau, args) * dtau; 
                }
                else if(bounds == ninf){
                     fx = f(b - tau, args) * dtau; 
                }
            }
            else{ // else proper integral
                fx = f(x, args); 
            }

            kronrodIntegral += fx * w; 
        }

        Complex gaussIntegral = 0; 
        for(int i = 0; i < gaussNodes.size(); i++){
            Complex x = ( (gaussNodes[i] + 1) / 2 ) * (bf - af) + af; 
            Complex w = gaussWeights[i] * ( (bf - af) / 2 );

            Complex fx; 
            if(bounds != proper) { // improper integral!
                Complex tau = (1 - x) / x; // substitution for change of integral bounds.. 
                Complex dtau = 1/(x*x); 
                if(bounds == twoSidedInf){
                    fx = (f(tau, args) + f(-tau, args)) * dtau; 
                }
                else if(bounds == pinf){
                    fx = f(a + tau, args) * dtau; 
                }
                else if(bounds == ninf){
                     fx = f(b - tau, args) * dtau; 
                }
            }
            else{ // else proper integral
                fx = f(x, args); 
            }

            gaussIntegral += fx * w; 
        }

        double epsilon = mod(kronrodIntegral - gaussIntegral); 

        if(epsilon > tau){
            Complex m = (af + bf) / 2; 
            kronrodIntegral = integral(f, af, m, tau/2, args) + integral(f, m, bf, tau/2, args); 
        }
        return kronrodIntegral;
    }





    // // ADD NEGATIVE INF. 
    // // AND ADD DNE LIMITS. 
    // constexpr Complex limit(Complex (*f)(Complex), const Complex& x) noexcept {
    //     if(x == INF) return f(1e30); 
    //     const double h = 1e-8;
    //     return f(x - h);
    // }

    // // This is a function that should only be used by the ratio test.
    // constexpr Complex limit(Complex (*ratioFn)(Complex, Complex (*)(Complex)), Complex (*f)(Complex), const Complex& x) noexcept {
    //     if(x == INF) return ratioFn(1e30, f); 
    //     const double h = 1e-8;
    //     return ratioFn(x - h, f);
    // }
}; 

#endif // NUMERICALANALYSIS_HPP