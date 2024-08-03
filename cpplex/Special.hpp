#ifndef SPECIAL_HPP
#define SPECIAL_HPP

#include <vector>
#include <utility>
#include "Complex.hpp"
#include "NumericalAnalysis.hpp"

namespace cpplex{

    class Special{


    }; 

    constexpr Complex gamma(const Complex& z) noexcept {
        Complex zf(z); 
        const int g = 8;
        const std::vector<double> C = {    
            0.9999999999999999298,
            1975.3739023578852322,
            -4397.3823927922428918,
            3462.6328459862717019,
            -1156.9851431631167820,
            154.53815050252775060,
            -6.2536716123689161798,
            0.034642762454736807441,
            -7.4776171974442977377e-7,
            6.3041253821852264261e-8,
            -2.7405717035683877489e-8,
            4.0486948817567609101e-9
        };

        if(zf.real() < 0.5){
            return M_PI / (sin(M_PI * zf) / gamma(1 - zf)); 
        }
        else{
            zf -= 1; // Lanczos is defined for gamma(z - 1)
            Complex A = 0; 
            for(int k = 1; k < C.size(); k++){
                A += C[k] / (zf + k);
            }
            A += C[0]; 

            return std::sqrt(2 * M_PI) * pow(zf + g + 0.5, zf + 0.5) * exp(-(zf +  g + 0.5)) * A; 
        }
    } 

    constexpr Complex loggamma(const Complex& z) noexcept {
        return log(gamma(z)); 
    }

    constexpr Complex rgamma(const Complex& z) noexcept {
        return 1 / gamma(z); 
    }

    // these all kinda suck for complex inputs
    // ALSO: ADD REAL INTEGRAL TOO SO I DONT HAVE TO KEEP DOING .real()
    constexpr double gammainc(const double a, const double x) noexcept {
        auto f = [](Complex t, std::vector<Complex> args) { return pow(t, args[0] - 1)  * exp(-t); };
        return integral(f, 0, x, {a}).real(); 
    }

    constexpr double gammaincc(const double a, const double x) noexcept {
        auto f = [](Complex t, std::vector<Complex> args) { return pow(t, args[0] - 1)  * exp(-t); };
        return integral(f, x, INF, {a}).real(); 
    }

    constexpr double gammaincReg(const double a, const double x) noexcept {
        return gammainc(a, x) / gamma(a).real();
    }

    constexpr double gammainccReg(const double a, const double x) noexcept {
        return gammaincc(a, x) / gamma(a).real();
    }

    constexpr double alpha(const int n, const double x) noexcept {
        auto f = [](Complex t, std::vector<Complex> args) { return pow(t, args[0])  * exp(-args[1] * t); };
        return integral(f, 1, INF, {n, x}).real(); 
    }

    constexpr double beta(const double a, const double b) noexcept {
        return (gamma(a) * gamma(b)).real() / gamma(a + b).real();
    }

    constexpr double betainc(const double x, const double a, const double b) noexcept {
        auto f = [](Complex t, std::vector<Complex> args) { return pow(t, args[0] - 1) * pow(1 - t, args[1] - 1); };
        return integral(f, 0, x, {a, b}).real();
    }

    constexpr double betaincReg(const double x, const double a, const double b) noexcept {
        return betainc(x, a, b) / beta(a, b); 
    }

    // rest of these are cool though 

    constexpr Complex pi(const Complex& z) noexcept {
        return gamma(z + 1); 
    }

    constexpr Complex psi(const Complex& z) noexcept {
        auto f = [](Complex z) { return gamma(z); };
        return derivative(f, z) / gamma(z);
    }

    constexpr Complex digamma(const Complex& z) noexcept {
        return psi(z); 
    }

    constexpr double zeta(const double x) noexcept {
        auto f = [](Complex t, std::vector<Complex> args) { return pow(t, args[0] - 1)  / (exp(t) - 1); };
        return (integral(f, 0, INF, {x}) / gamma(x)).real(); 
    }

    // This is only for integer order n
    constexpr Complex jn(const int n, const Complex& z) noexcept {
        auto f = [](Complex theta, std::vector<Complex> args) { return cos(args[0] * sin(theta) - args[1] * theta); };
        return (1.0 / M_PI) * integral(f, 0, M_PI, {z, n});
    }

    // https://dlmf.nist.gov/10.9
    inline Complex jv(const double v, const Complex& z) noexcept {
        auto f = [](Complex theta, std::vector<Complex> args) { return cos(args[0] * cos(theta)) * pow(sin(theta), 2 * args[1]); };
        Complex multiplier = pow(0.5 * z, v) / (std::sqrt(M_PI) * gamma(v + 0.5)); 
        return multiplier * integral(f, 0, M_PI, {z, v});
    }
    // https://dlmf.nist.gov/10.9
    inline Complex yv(const double v, const Complex& z) noexcept {
        auto f1 = [](Complex theta, std::vector<Complex> args) { return pow(1 - theta * theta, args[1] - 0.5) * sin(args[0] * theta); };
        auto f2 = [](Complex theta, std::vector<Complex> args) { return exp(-args[0] * theta) * pow(1 + theta * theta, args[1] - 0.5); };
        Complex multiplier = 2 * pow(0.5 * z, v) / (std::sqrt(M_PI) * gamma(v + 0.5)); 
        return multiplier * ( integral(f1, 0, 1, {z, v}) - integral(f2, 0, INF, {z, v}) ); 
    }

    constexpr Complex hv1(const double v, const Complex& z) noexcept {
        return jv(v, z) + 1_j * yv(v, z); 
    }

    constexpr Complex hv2(const double v, const Complex& z) noexcept {
        return jv(v, z) - 1_j * yv(v, z); 
    }

    // constexpr Complex hyp2f1(const double a, const double b, const double c, const double z) noexcept {
    //     auto f = [](Complex t, std::vector<Complex> args) {return pow(t, args[1] - 1) * pow(1 - t, args[2] - args[1] - 1) / pow(1 - args[3] * t, -args[0]); };
    //     return (1.0 / beta(b, c - b)) * integral(f, 0, 1, {a, b, c, z});
    // }

    // change to constexpr later. 
    // probably will just switch this out for cpplex::sqrt... but it would be more expensive, technically. 
    inline Complex erf(const Complex& z) noexcept {

        auto f = [](Complex z) { return exp(-z*z); };

        return (2/std::sqrt(M_PI)) * integral(f, 0, z); 

        // double x = z.getX(); 
        // double y = z.getY(); 

        // auto f1 = [](Complex z, std::vector<Complex> args) { return exp(z*z) * cos(2 * z * args[0]); };
        // auto f2 = [](Complex z, std::vector<Complex> args) { return exp(z*z) * sin(2 * z * args[0]); };

        // double norm = 2 / std::sqrt(M_PI) * std::exp(-x * x);
        // return std::erf(x) + norm * Complex(integral(f2, 0, y, {x}).real(), integral(f1, 0, y, {x}).real()); 
    }

    constexpr Complex erfc(const Complex& z) noexcept {
        return 1 - erf(z); 
    }

    constexpr Complex erfcx(const Complex& z) noexcept {
        return exp(z * z) * erfc(z); 
    }

    // change to constexpr later. 
    inline Complex erfi(const Complex& z) noexcept {

        auto f = [](Complex z) { return exp(z*z); };

        return (2/std::sqrt(M_PI)) * integral(f, 0, z); 
    }

    constexpr Complex wofz(const Complex& z) noexcept {
        return exp(-z * z) * erfc(-1_j * z); 
    }

    constexpr Complex dawsn(const Complex& z) noexcept {
        auto f = [](Complex z) { return exp(z*z); };
        return exp(-z * z) * integral(f, 0, z);
    }

    constexpr std::pair<Complex, Complex> fresnel(const Complex& z) noexcept {
        auto s = [](Complex z) { return sin(M_PI * z*z / 2); };
        auto c = [](Complex z) { return cos(M_PI * z*z / 2); };
        return {integral(s, 0, z), integral(c, 0, z)};
    }

    constexpr Complex li(const Complex& z) noexcept {

        auto f = [](Complex z) { return 1 / log(z); };

        return integral(f, 2, z) + li2; 
    }

    constexpr Complex ei(const Complex& z) noexcept {

        auto f = [](Complex z) { return (exp(z) - 1) / z; };

        return integral(f, 0, z) + 0.5 * (log(z) - log(1 / z)) + eulerGamma; 
    }

    constexpr Complex productlog(const Complex& z, double epochs=1000) noexcept {
        Complex w = 0; 

        for(int i = 0; i < epochs; i++){
            w -= (w * exp(w) - z) / (exp(w) + w * exp(w)); 
        }

        return w; 
    }

    constexpr Complex weightomega(const Complex& z, double epochs=1000) noexcept {
        Complex w = z; // just a guess 

        for(int i = 0; i < epochs; i++){
            w -= (w + log(w) - z) / (1 + 1/w); 
        }

        return w; 
    }

    // Cheyb. Polys, 1 + 2. 
    constexpr Complex tn(const int n, Complex z) noexcept {
        return cos(n * acos(z)); 
    }

    constexpr Complex un(const int n, Complex z) noexcept {
        return (sin((n + 1) * acos(z))) / (sin(acos(z)));
    }

    constexpr Complex sigmoid(const Complex& z) noexcept {
        return 1.0 / (1 + exp(-z)); 
    }

    constexpr double sgn(const double x){
        if(x < 0) return -1; 
        else if(x == 0) return 0; 
        else return 1; 
    }

    constexpr int factorial(const int k){
        int fact = 1; 
        for(int i = 1; i <= k; i++) fact *= i; 
        return fact; 
    }

    constexpr int binom(const int n, const int k){
        return factorial(n) / (factorial(k) * factorial(n - k)); 
    }
}

#endif // SPECIAL_HPP