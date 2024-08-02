#ifndef FFT_HPP
#define FFT_HPP

#include <vector> 
#include <utility> 
#include <cmath>
#include "Complex/Complex.hpp"

namespace cpplex {

     // constexpr std::vector<Complex> fft(const std::vector<Complex>& X) noexcept; 
    class FFT {



    }; 

    constexpr std::pair<std::vector<double>, std::vector<double>> decomp(const std::vector<Complex>& X) noexcept {   
        int n = X.size(); 
        std::vector<double> reX(n); 
        std::vector<double> imX(n); 

        for(int i = 0; i < n; i++){
            reX[i] = X[i].getX();
            imX[i] = X[i].getY();
        }
        return {reX, imX}; 
    }

    constexpr std::vector<Complex> conj(const std::vector<Complex>& X) noexcept {   
        int n = X.size(); 
        std::vector<Complex> conjX(n);
        for(int i = 0; i < n; i++){
            conjX[i] = conj(X[i]); 
        }
        return conjX;
    }

    inline namespace{ 

        constexpr std::vector<Complex> fft_(const std::vector<Complex>& X) noexcept {       
            int n = X.size(); 

            if(n == 1) return X; 

            Complex omega = exp(2 * M_PI * 1_j / double(n)); 
            std::vector<Complex> E(n / 2); 
            std::vector<Complex> O(n / 2); 

            for(int i = 0; i < n/2; i++){
                E[i] = X[2 * i]; 
                O[i] = X[2 * i + 1]; 
            }

            std::vector<Complex> Ye = fft_(E);
            std::vector<Complex> Yo = fft_(O);

            std::vector<Complex> Y(n); 

            for(int i = 0; i < n / 2; i++){
                Complex powOmega = pow(omega, i); 
                Y[i] = Ye[i] + powOmega * Yo[i]; 
                Y[i + n/2] = Ye[i] - powOmega * Yo[i]; 
            }
            return Y; 

        }

        constexpr std::vector<Complex> ifft_(const std::vector<Complex>& X) noexcept {
            int n = X.size(); 
            
            std::vector<Complex> Y = fft_(conj(X)); 
            for(int i = 0; i < n; i++){
                Y[i] = conj(Y[i]) / n; 
            }
            return Y; 

        }

        constexpr std::vector<Complex> czt_(const std::vector<Complex>& X) noexcept{
            double n = X.size();
            double m = std::pow(2, std::ceil(std::log2(2 * n - 1))); // next power of two. 
            Complex omega = exp(-M_PI * 1_j / double(n)); 

            std::vector<Complex> a(m);
            for(double i = 0; i < n; i++){
                a[i] = pow(omega, i*i) * X[i]; 
            }

            std::vector<Complex> b(m);
            for(double i = 0; i < n; i++){
                b[i] = pow(omega, -i*i);
                if (i > 0) b[m-i] = b[i];
            }


            std::vector<Complex> IFFT(m); 
            std::vector<Complex> aFFT = fft_(a); 
            std::vector<Complex> bFFT = fft_(b); 

            for(int i = 0; i < m; i++){
                IFFT[i] = aFFT[i] * bFFT[i]; 
            }

            IFFT = ifft_(IFFT); 

            std::vector<Complex> Y(n);
            for(int i = 0; i < n; i++){
                Y[i] = conj(b[i]) * IFFT[i];
            }
            return Y; 
        }
    }

    constexpr std::vector<Complex> fft(const std::vector<Complex>& X) noexcept {       
        int n = X.size(); 

        double log2n = std::log2(n); 
        if(int(floor(log2n)) != log2n) {
            return czt_(X); 
        }
        else{
            std::vector<Complex> Y = fft_(X); 

            std::vector<Complex> Yf(Y);

            for(int i = 1; i < n; i++){
                Yf[i] = Y[n - i];
            }
            return Yf; 
        }
    }


    constexpr std::vector<Complex> ifft(const std::vector<Complex>& X) noexcept {       
        int n = X.size(); 

        std::vector<Complex> conjX(X);
        for(int i = 0; i < n; i++){
            conjX[i] = conj(X[i]); 
        }
        std::vector<Complex> Y = fft(conjX); 
        for(int i = 0; i < n; i++){
            Y[i] = conj(Y[i]) / n; 
        }
        return Y; 
    }

    inline namespace {
        constexpr std::vector<Complex> dct2_(const std::vector<Complex>& X) noexcept {
                int n = X.size(); 
                std::vector<Complex> Y(2 * n); 
                for(int i = 0; i < n; i++){
                    Y[i] = X[i];
                }
                for(int i = 0; i < n; i++){
                    Y[i + n] = X[(n - 1) - i]; 
                }

                Y = fft(Y); 

                std::vector<Complex> Yf(n);

                for(int i = 0; i < n; i++){
                    Yf[i] = Y[i] * exp(-1_j * M_PI * double(i) / double(2 * n));
                    Yf[i] = Yf[i].real();

                }

                return Yf; 
        }

        constexpr std::vector<Complex> dct3_(const std::vector<Complex>& X) noexcept {

            int n = X.size(); 
            std::vector<Complex> Y(n); 

            // Multiply by 2 * n because DCT 3 is not the inverse of DCT2; in fact, it is the inverse of the DCT 2 * (2n)
            for(int i = 0; i < n; i++){
                Y[i] = X[i] * exp(1_j * M_PI * double(i) / double(2 * n)) * (2 * n); 
            }
            Y[0] /= 2.0;  // as per formula

            Y = ifft(Y); 

            for(int i = 0; i < n; i++){
                Y[i] = Y[i].real(); 
            }

            std::vector<Complex> Yf(n); 
            for(int i = 0; i < n/2; i++){
                Yf[2 * i] = Y[i]; 
                Yf[2 * i + 1] = Y[n - 1 - i]; 
            }
            if(n % 2 != 0) Yf[n - 1] = Y[n/2]; 

            return Yf; 
        }


        // https://github.com/ARM-software/CMSIS_4/blob/master/CMSIS/DSP_Lib/Source/TransformFunctions/arm_dct4_f32.c
        // See comments explaining the algorithm lines 149 --> 166. 
        constexpr std::vector<Complex> dct4_(const std::vector<Complex>& X) noexcept {

            int n = X.size(); 
            std::vector<Complex> Y(n); 

            for(int i = 0; i < n; i++){
                Y[i] = 2 * X[i] * std::cos(double(M_PI * (2 * i + 1)) / double(4 * n));
            }

            Y = dct2_(Y); // we already have a working dct2, so no need to re-implement using FFT. 

            std::vector<Complex> Yf(n);
            Yf[0]  = 0.5 * Y[0]; 
            for(int i = 1; i < n; i++){
                Yf[i] = Y[i] - Yf[i - 1]; 
            }

            return Yf; 
        }
    }

    // TO DO: Throw an error if wrong type. 

    constexpr std::vector<Complex> dct(const std::vector<Complex>& X, int type = 2) noexcept {
        if(type == 1){
            int n = X.size();
            std::vector<Complex> reX(2 * n - 2); 
            std::vector<Complex> imX(2 * n - 2); 

            for(int i = 0; i < n; i++){
                reX[i] = X[i].getX();
                imX[i] = X[i].getY();
            }
            for(int i = n; i < 2 * n - 2; i++){
                reX[i] = X[(n - 2) - (i - n)].getX(); 
                imX[i] = X[(n - 2) - (i - n)].getY(); 
            }

            std::vector<Complex> reY = fft(reX); 
            std::vector<Complex> imY = fft(imX); 

            std::vector<Complex> Y(n); 
            for(int i = 0; i < n; i++){
                Y[i] = reY[i].getX() + imY[i].getX() * 1_j; 
            }

            return Y; 
        }
        else if(type == 2){

            int n = X.size();
            std::vector<Complex> reX(n); 
            std::vector<Complex> imX(n); 

            for(int i = 0; i < n; i++){
                reX[i] = X[i].getX();
                imX[i] = X[i].getY();
            }

            std::vector<Complex> reY = dct2_(reX); 
            std::vector<Complex> imY = dct2_(imX); 

            std::vector<Complex> Y(n); 
            for(int i = 0; i < n; i++){
                Y[i] = reY[i] + imY[i] * 1_j; 
            }

            return Y; 

        }
        else if(type == 3){
            int n = X.size();
            std::vector<Complex> reX(n); 
            std::vector<Complex> imX(n); 

            for(int i = 0; i < n; i++){
                reX[i] = X[i].getX();
                imX[i] = X[i].getY();
            }

            std::vector<Complex> reY = dct3_(reX); 
            std::vector<Complex> imY = dct3_(imX); 

            std::vector<Complex> Y(n); 
            for(int i = 0; i < n; i++){
                Y[i] = reY[i] + imY[i] * 1_j; 
            }

            return Y; 

            
        }
        // https://www.appletonaudio.com/blog/2013/derivation-of-fast-dct-4-algorithm-based-on-dft/
        else if(type == 4){
            int n = X.size();
            std::vector<Complex> reX(n); 
            std::vector<Complex> imX(n); 

            for(int i = 0; i < n; i++){
                reX[i] = X[i].getX();
                imX[i] = X[i].getY();
            }

            std::vector<Complex> reY = dct4_(reX); 
            std::vector<Complex> imY = dct4_(imX); 

            std::vector<Complex> Y(n); 
            for(int i = 0; i < n; i++){
                Y[i] = reY[i] + imY[i] * 1_j; 
            }

            return Y; 
        }

    }

    constexpr std::vector<Complex> idct(const std::vector<Complex>& X, int type=2) noexcept {
        if(type == 1) return dct(X, 1); 
        else if(type == 2) return dct(X, 3); 
        else if(type == 3) return dct(X, 2); 
        else if(type == 4) return dct(X, 4); 
    }

    inline namespace {

        inline std::vector<Complex> dst2_(const std::vector<Complex>& X) noexcept {
            int n = X.size(); 
            std::vector<Complex> negX(X); 

            for(int i = 0; i < n; i++){
                negX[i] *= std::pow(-1, i);
            }

            std::vector<Complex> Y = dct(negX, 2); 


            std::vector<Complex> Yf(n); 
            for(int i = 0; i < n; i++){
                Yf[i] = Y[n - 1 - i]; 
            }
            return Yf;
        }    

        inline std::vector<Complex> dst3_(const std::vector<Complex>& X) noexcept {
            int n = X.size(); 
            std::vector<Complex> revX(n); 

            for(int i = 0; i < n; i++){
                revX[i] = X[n - 1 - i]; 
            }

            std::vector<Complex> Y = dct(revX, 3); 

            for(int i = 0; i < n; i++){
                Y[i] *= std::pow(-1, i); 
            }
            return Y;
        }    

        inline std::vector<Complex> dst4_(const std::vector<Complex>& X) noexcept {
            int n = X.size(); 
            std::vector<Complex> negX(X); 

            for(int i = 0; i < n; i++){
                negX[i] *= std::pow(-1, i);
            }

            std::vector<Complex> Y = dct(negX, 4); 


            std::vector<Complex> Yf(n); 
            for(int i = 0; i < n; i++){
                Yf[i] = Y[n - 1 - i]; 
            }
            return Yf;
        }    
    }

    constexpr std::vector<Complex> dst(const std::vector<Complex>& X, int type = 2) noexcept {
        if(type == 1){
            int n = X.size();
            std::vector<Complex> reX(2 * n + 2); 
            std::vector<Complex> imX(2 * n + 2); 

            for(int i = 1; i < n + 1; i++){
                reX[i] = X[i - 1].getX();
                imX[i] = X[i - 1].getY();
            }
            for(int i = n + 2; i < 2 * n + 2; i++){
                reX[i] = -X[2 * n + 1 - i].getX(); 
                imX[i] = -X[2 * n + 1 - i].getY(); 
            }

            std::vector<Complex> reY = fft(reX); 
            std::vector<Complex> imY = fft(imX); 

            std::vector<Complex> Y(n); 
            for(int i = 2 * n + 1; i >= n + 2; i--){
                Y[2 * n + 1 - i] = reY[i].getY() + imY[i].getY() * 1_j; 
            }

            return Y; 
        }
        else if(type == 2){

            int n = X.size();
            std::vector<Complex> reX(n); 
            std::vector<Complex> imX(n); 

            for(int i = 0; i < n; i++){
                reX[i] = X[i].getX();
                imX[i] = X[i].getY();
            }

            std::vector<Complex> reY = dst2_(reX); 
            std::vector<Complex> imY = dst2_(imX); 

            std::vector<Complex> Y(n); 
            for(int i = 0; i < n; i++){
                Y[i] = reY[i] + imY[i] * 1_j; 
            }

            return Y; 

        }
        else if(type == 3){
            int n = X.size();
            std::vector<Complex> reX(n); 
            std::vector<Complex> imX(n); 

            for(int i = 0; i < n; i++){
                reX[i] = X[i].getX();
                imX[i] = X[i].getY();
            }

            std::vector<Complex> reY = dst3_(reX); 
            std::vector<Complex> imY = dst3_(imX); 

            std::vector<Complex> Y(n); 
            for(int i = 0; i < n; i++){
                Y[i] = reY[i] + imY[i] * 1_j; 
            }

            return Y; 

            
        }
        // https://www.appletonaudio.com/blog/2013/derivation-of-fast-dct-4-algorithm-based-on-dft/
        else if(type == 4){
            int n = X.size();
            std::vector<Complex> reX(n); 
            std::vector<Complex> imX(n); 

            for(int i = 0; i < n; i++){
                reX[i] = X[i].getX();
                imX[i] = X[i].getY();
            }

            std::vector<Complex> reY = dst4_(reX); 
            std::vector<Complex> imY = dst4_(imX); 

            std::vector<Complex> Y(n); 
            for(int i = 0; i < n; i++){
                Y[i] = reY[i] + imY[i] * 1_j; 
            }

            return Y; 
        }

    }

    constexpr std::vector<Complex> idst(const std::vector<Complex>& X, int type=2) noexcept {
        if(type == 1) return dst(X, 1); 
        else if(type == 2) return dst(X, 3); 
        else if(type == 3) return dst(X, 2); 
        else if(type == 4) return dst(X, 4); 
    }

    inline std::vector<Complex> dcht(const std::vector<Complex>& X) noexcept {
        int n = X.size(); 
        std::vector<Complex> Y = dct(X); 
        for(int i = 0; i < n; i++){
            Y[i] *= std::sqrt(2.0 / n);
        }
        Y[0] /= std::sqrt(2);
        return Y;
    }

    inline std::vector<Complex> idcht(const std::vector<Complex>& X) noexcept {
        int n = X.size(); 
        std::vector<Complex> Y(X);
        for(int i = 0; i < n; i++){
            Y[i] *= std::sqrt(n / 2.0);
        }
        Y[0] /= std::sqrt(2);
        return idct(Y); 
    }

    // discrete hartley transform.. 
    constexpr std::vector<Complex> dhat(const std::vector<Complex>& X) noexcept {
            int n = X.size();
            std::vector<Complex> reX(n); 
            std::vector<Complex> imX(n); 

            for(int i = 0; i < n; i++){
                reX[i] = X[i].getX();
                imX[i] = X[i].getY();
            }

            std::vector<Complex> reY(n); 
            std::vector<Complex> imY(n); 

            reY = fft(reX);
            imY = fft(imX);

            std::vector<Complex> Y(n);
            for(int i = 0; i < n; i++){
                Y[i] = (reY[i] * (1 + 1_j)).real() + (imY[i] * (1 + 1_j)).real() * 1_j; 
            }

            return Y;  

    }

    constexpr std::vector<Complex> idhat(const std::vector<Complex>& X) noexcept {
            int n = X.size();
            std::vector<Complex> Y = dhat(X);

            for(int i = 0; i < n; i++){
                Y[i] /= n; 
            }
            return Y;
    }

    inline namespace {

        constexpr std::vector<Complex> dht_(const std::vector<Complex>& X) noexcept {
            int n = X.size(); 
            std::vector<Complex> Y = fft(X); 
            std::vector<Complex> H(n); 

            if(n % 2 == 0){ // even
                for(int i = 0; i < n; i++){
                    if(i == 0){
                        H[i] = 1; 
                    }
                    else if(i >= 1 && i <= n/2 - 1){
                        H[i] = 2; 
                    }
                    else if(i == n / 2){
                        H[i] = 1; 
                    }
                    else{
                        H[i] = 0; 
                    }
                }
            }
            else{ // odd 
                for(int i = 0; i < n; i++){
                    if(i == 0){
                        H[i] = 1; 
                    }
                    else if(i >= 1 && i <= (n-1)/2){
                        H[i] = 2; 
                    }
                    else{
                        H[i] = 0; 
                    }
                }
            }  

            for(int i = 0; i < n; i++){
                Y[i] *= H[i]; 
            }
            return ifft(Y);

        }

        constexpr std::vector<Complex> idht_(const std::vector<Complex>& X) noexcept {
            int n = X.size(); 
            std::vector<Complex> Y = fft(X); 
            std::vector<Complex> H(n); 

            if(n % 2 == 0){ // even
                for(int i = 0; i < n; i++){
                    if(i == 0){
                        H[i] = 1; 
                    }
                    else if(i >= 1 && i <= n/2 - 1){
                        H[i] = 0.5; 
                    }
                    else if(i == n / 2){
                        H[i] = 1; 
                    }
                    else{
                        H[i] = 0; 
                    }
                }
            }
            else{ // odd 
                for(int i = 0; i < n; i++){
                    if(i == 0){
                        H[i] = 1; 
                    }
                    else if(i >= 1 && i <= (n-1)/2){
                        H[i] = 0.5; 
                    }
                    else{
                        H[i] = 0; 
                    }
                }
            }  

            for(int i = 0; i < n; i++){
                Y[i] *= H[i]; 
            }
            return ifft(Y);
        }

    }

    // hilbert transform. 
    // https://dsp.stackexchange.com/questions/29263/is-hilbert-transform-not-defined-for-complex-signals
    constexpr std::vector<Complex> dht(const std::vector<Complex>& X) noexcept {
        int n = X.size();
        std::vector<Complex> reX(n); 
        std::vector<Complex> imX(n); 

        bool real = true; 
        for(int i = 0; i < n; i++){
            reX[i] = X[i].getX();
            imX[i] = X[i].getY();
            if(imX[i] != 0) real = false; 
        }
        if(real) return dht_(X);

        std::vector<Complex> reY = dht_(reX); 
        std::vector<Complex> imY = dht_(imX); 
        std::vector<Complex> Yf(n); 

        for(int i = 0; i < n; i++){
            Yf[i] = reY[i].im() + 1_j * imY[i].im(); 
        }
        return Yf;
    }

    constexpr std::vector<Complex> idht(const std::vector<Complex>& X) noexcept {
        int n = X.size();
        std::vector<Complex> reX(n); 
        std::vector<Complex> imX(n); 

        bool real = true; 
        for(int i = 0; i < n; i++){
            reX[i] = X[i].getX();
            imX[i] = X[i].getY();
            if(imX[i] != 0) real = false; 
        }
        if(real) return idht_(X);

        std::vector<Complex> reY = idht_(reX); 
        std::vector<Complex> imY = idht_(imX); 
        std::vector<Complex> Yf(n); 

        for(int i = 0; i < n; i++){
            Yf[i] = reY[i].im() + 1_j * imY[i].im(); 
        }
        return Yf;
    }

//     // we need a separte fn because czt_ needs to use fft_ and ifft_ which only handle power of 2 input arrays. 
//     constexpr std::vector<Complex> zt(const std::vector<Complex>& X, const Complex& omega, const double m, const double alpha) noexcept {
//         double n = X.size();

//         std::vector<Complex> a(m);
//         for(double i = 0; i < n; i++){
//             a[i] = pow(omega * alpha, i*i) * X[i]; 
//         }

//         std::vector<Complex> b(m);
//         for(double i = 0; i < n; i++){
//             b[i] = pow(omega * alpha, -i*i);
//             if (i > 0) b[m-i] = b[i];
//         }


//         std::vector<Complex> IFFT(m); 
//         std::vector<Complex> aFFT = fft(a); 
//         std::vector<Complex> bFFT = fft(b); 

//         for(int i = 0; i < m; i++){
//             IFFT[i] = aFFT[i] * bFFT[i]; 
//         }

//         IFFT = ifft(IFFT); 

//         std::vector<Complex> Y(m);
//         for(int i = 0; i < m; i++){
//             Y[i] = conj(b[i]) * IFFT[i];
//         }
//         return Y; 
//     }
}

#endif // FFT_HPP