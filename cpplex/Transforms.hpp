#ifndef TRANSFORMS_HPP
#define TRANSFORMS_HPP

#include <vector> 
#include <cmath>
#include "Complex.hpp"
#include "Special.hpp"

namespace cpplex {
    constexpr Complex laplace(Complex (*f)(Complex), const Complex& z, const Complex& a = 1) noexcept {
        auto ft = [](Complex t, std::vector<Complex (*)(Complex)> fargs, std::vector<Complex> args) { return fargs[0](args[0] * t) * exp(-t * args[1]); };
        return integral(ft, 0, INF, {f}, {a, z}); 
    }

    // these all do not work because they do not converge.
    // constexpr Complex hankel(Complex (*f)(Complex), const double v, const Complex& z){
    //     auto ft = [](Complex t, std::vector<Complex (*)(Complex)> fargs, std::vector<Complex> args) { return fargs[0](t) * jv(args[0].real(), args[1] * t) * t; };
    //     return integral(ft, 0, INF, {f}, {v, z});       
    // }

    // constexpr Complex fourier(Complex (*f)(Complex), const Complex& z) noexcept {
    //     auto ft = [](Complex t, Complex (*f)(Complex), std::vector<Complex> args) { return f(t) * exp(-1_j * 2 * M_PI * t * args[0]); };
    //     return integral(ft, -INF, INF, f, {z}); 
    // }

    // constexpr Complex invfourier(Complex (*f)(Complex), const Complex& z) noexcept {
    //     auto ft = [](Complex t, Complex (*f)(Complex), std::vector<Complex> args) { return f(t) * exp(1_j * 2 * M_PI * t * args[0]); };
    //     return integral(ft, -INF, INF, f, {z}); 
    // }

    // inline Complex mellin(Complex (*f)(Complex), const Complex& z) noexcept {
    //     auto ft = [](Complex t, Complex (*f)(Complex), std::vector<Complex> args) {return f(t) * pow(t, args[0] - 1);  };
    //     return integral(ft, 0, INF, f, {z});
    // 
    // }

    inline Complex hartley(Complex (*f)(Complex), const Complex& z) noexcept {
        auto ft = [](Complex t, std::vector<Complex (*)(Complex)> fargs, std::vector<Complex> args) { return fargs[0](t) * (sin(args[0] * t) + cos(args[0] * t));  };
        return integral(ft, -INF, INF, {f}, {z}) / std::sqrt(2 * M_PI); 
    }

    // Hartley is its own inverse. W. 
    constexpr Complex invhartley(Complex (*f)(Complex), const Complex& z) noexcept {
        return hartley(f, z); 
    }

    constexpr Complex conv(Complex (*f)(Complex), Complex (*g)(Complex), const Complex& z) noexcept {
        auto ft = [](Complex t, std::vector<Complex (*)(Complex)> fargs, std::vector<Complex> args) { return fargs[0](args[0] - t) * fargs[1](t);  };
        return integral(ft, NINF, INF, {f, g}, {z}); 
    }

    constexpr Complex crosscorr(Complex (*f)(Complex), Complex (*g)(Complex), const Complex& z) noexcept {
        auto ft = [](Complex t, std::vector<Complex (*)(Complex)> fargs,  std::vector<Complex> args) { 
            Complex fz = fargs[0](t); 
            return (fz.real() - fz.im() * 1_j) * fargs[1](t + args[0]); 
        };
        return integral(ft, NINF, INF, {f, g}, {z}); 

    }

    constexpr Complex autocorr(Complex (*f)(Complex), const Complex& z) noexcept {
        return crosscorr(f, f, z); 
    }   

    constexpr Complex inner(Complex (*f)(Complex), Complex (*g)(Complex), const Complex& z, const Complex& a, const Complex& b) noexcept {
        auto ft = [](Complex t, std::vector<Complex (*)(Complex)> fargs) { 
            Complex gz = fargs[1](t); 
            return fargs[0](t) * (gz.real() - gz.im() * 1_j); 
        };
        return integral(ft, a, b, {f, g}); 

    }

}


#endif // TRANSFORMS_HPP