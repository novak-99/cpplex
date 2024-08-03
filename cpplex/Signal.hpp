#ifndef SIGNAL_HPP
#define SIGNAL_HPP

#include <vector> 
#include <cmath>
#include <algorithm> 
#include "Complex.hpp"
#include "Transforms.hpp"

namespace cpplex {

    // Thank you! https://dsp.stackexchange.com/questions/736/how-do-i-implement-cross-correlation-to-prove-two-audio-files-are-similar

    constexpr std::vector<Complex> conv(const std::vector<Complex>& X, const std::vector<Complex>& Y) noexcept {       
        int n = X.size(); 
        int m = Y.size(); 
        int l = n + m - 1; 
        std::vector<Complex> padX(l); 
        std::vector<Complex> padY(l); 

        std::copy(X.begin(), X.begin() + n, padX.begin());
        std::copy(Y.begin(), Y.begin() + m, padY.begin());

        std::vector<Complex> fftX = fft(padX);
        std::vector<Complex> fftY = fft(padY);

        std::vector<Complex> C(l); 

        for(int i = 0; i < l; i++){
            C[i] = fftX[i] * fftY[i]; 
        }

        return ifft(C); 
    }

    constexpr std::vector<Complex> crosscorr(const std::vector<Complex>& X, const std::vector<Complex>& Y) noexcept {       
        int n = X.size(); 
        int m = Y.size(); 
        int l = n + m - 1; 
        std::vector<Complex> padX(l); 
        std::vector<Complex> padY(l); 

        std::vector<Complex> revY(m); 
        for(int i = 0; i < m; i++){
            revY[i] = Y[m - 1 - i];
        }

        std::copy(X.begin(), X.begin() + n, padX.begin());
        std::copy(revY.begin(), revY.begin() + m, padY.begin());

        std::vector<Complex> fftX = fft(padX);
        std::vector<Complex> fftConjY = fft(conj(padY));

        std::vector<Complex> C(l); 

        for(int i = 0; i < l; i++){
            C[i] = fftX[i] * fftConjY[i]; 
        }

        return ifft(C); 
    }

    constexpr std::vector<Complex> autocorr(const std::vector<Complex>& X) noexcept {       
        return crosscorr(X, X); 
    }

    constexpr Complex inner(const std::vector<Complex>& X, const std::vector<Complex>& Y) noexcept {
        Complex sum = 0; 
        int n = X.size(); 
        for(int i = 0; i < n; i++) sum += X[i] * conj(Y[i]);
        return sum; 
    }
}

#endif // SIGNAL_HPP