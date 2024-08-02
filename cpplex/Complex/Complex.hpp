#ifndef COMPLEX_HPP
#define COMPLEX_HPP

#include <iostream>
#include <complex> // only to convert from cpplex -> std, std -> cpplex. 
#include <vector> // for q roots of unity. 
#include <string>

namespace cpplex{

    class Complex{
        public:
            constexpr Complex() noexcept : x(0), y(0) {

            }

            constexpr Complex(const double x) noexcept : x(x), y(0) {

            }
            
            constexpr Complex(const double x, const double y) noexcept : x(x), y(y) {

            }

            constexpr Complex(const Complex& z) noexcept : x(z.getX()), y(z.getY()) {

            }

            template <class T>
            constexpr Complex(const std::complex<T>& z) noexcept : x(z.real()), y(z.imag()) {

            }

            constexpr void setX(const double x) noexcept {
                this->x = x; 
            }

            constexpr void setY(const double y) noexcept {
                this->y = y; 
            }

            constexpr double getX() const noexcept {
                return x; 
            }

            constexpr double getY() const noexcept {
                return y; 
            }

            constexpr double real() const noexcept {
                return x; 
            }

            constexpr double im() const noexcept {
                return y; 
            }

            constexpr Complex& operator=(const Complex& z) noexcept {
                x = z.getX(); 
                y = z.getY(); 
                return *this; 
            }

            constexpr Complex& operator=(const double alpha) noexcept {
                x = alpha; 
                y = 0; 
                return *this; 
            }

            constexpr Complex& operator+=(const Complex& z) noexcept {
                x += z.getX(); 
                y += z.getY(); 
                return *this; 
            }

            constexpr Complex& operator-=(const Complex& z) noexcept {
                x -= z.getX(); 
                y -= z.getY(); 
                return *this; 
            }

            constexpr Complex& operator*=(const Complex& z) noexcept {
                double a = x;
                double b = y;
                double c = z.getX();
                double d = z.getY(); 

                x = a * c - b * d; 
                y = a * d + c * b; 

                return *this; 
            }

            constexpr Complex& operator/=(const Complex& z) noexcept {
                double a = x; 
                double b = y;
                double c = z.getX();  
                double d = z.getY(); 

                double den = c * c + d * d; 
                x = (a * c + b * d)/den; 
                y = (b * c - a * d)/den; 


                // // Yeah, so this basically only causes a problem with atan2(x). 
                // // Remove this later & update atan2(x).
                // if(x == -0) { 
                //     x = 0; 
                // }
                // if(y == -0) { 
                //     y = 0; 
                // }
                return *this; 
            }

            constexpr Complex& operator+=(const double alpha) noexcept{
                x += alpha; 
                return *this;  
            }

            constexpr Complex& operator-=(const double alpha) noexcept{
                x -= alpha; 
                return *this;  
            }

            constexpr Complex& operator*=(const double alpha) noexcept{
                x *= alpha; 
                y *= alpha;
                return *this; 
            }

            constexpr Complex& operator/=(const double alpha) noexcept{
                x /= alpha; 
                y /= alpha;
                return *this; 
            }

            constexpr std::complex<double> std() noexcept {
                return {this->getX(), this->getY()};
            }

        private:
            double x; 
            double y; 
    };

    constexpr Complex operator+(const Complex& z) noexcept {
        return z; 
    }

    constexpr Complex operator-(const Complex& z) noexcept {
        return Complex(-z.getX(), -z.getY());
    }

    constexpr Complex operator+(const Complex& z1, const Complex& z2) noexcept { 
        Complex zf(z1);
        zf += z2;

        return zf;
    }

    constexpr Complex operator-(const Complex& z1, const Complex& z2) noexcept { 
        Complex zf(z1);
        zf -= z2;

        return zf;
    }

    constexpr Complex operator*(const Complex& z1, const Complex& z2) noexcept { 
        Complex zf(z1);
        zf *= z2;

        return zf;
    }

    constexpr Complex operator/(const Complex& z1, const Complex& z2) noexcept { 
        Complex zf(z1);
        zf /= z2;

        return zf;
    }

    constexpr Complex operator+(const Complex& z, const double alpha) noexcept {
        Complex zf(z);
        zf += alpha;

        return zf;
    }

    constexpr Complex operator-(const Complex& z, const double alpha) noexcept {
        Complex zf(z);
        zf -= alpha;

        return zf;
    }

    constexpr Complex operator*(const Complex& z, const double alpha) noexcept {
        Complex zf(z);
        zf *= alpha;

        return zf;
    }

    constexpr Complex operator/(const Complex& z, const double alpha) noexcept {
        Complex zf(z);
        zf /= alpha;

        return zf;
    }

    constexpr Complex operator+(const double alpha, const Complex& z) noexcept {
        Complex zf(z);
        zf += alpha;

        return zf;
    }

    constexpr Complex operator-(const double alpha, const Complex& z) noexcept {
        Complex zf(-z);
        zf += alpha;

        return zf;
    }

    constexpr Complex operator*(const double alpha, const Complex& z) noexcept {
        Complex zf(z);
        zf *= alpha;

        return zf;
    }

    constexpr Complex operator/(const double alpha, const Complex& z) noexcept {
        Complex zf(alpha);
        zf /= z; 

        return zf;
    }

    constexpr bool operator==(const Complex& z1, const Complex& z2) noexcept {
        return true ? (z1.getX() == z2.getX() && z1.getY() == z2.getY()) : false; 
    }

    constexpr bool operator==(const Complex& z, const double alpha) noexcept {
        return true ? (z.getX() == alpha && z.getY() == 0) : false; 
    }

    constexpr bool operator==(const double alpha, const Complex& z) noexcept {
        return true ? (z.getX() == alpha && z.getY() == 0) : false; 
    }

    constexpr bool operator!=(const Complex& z1, const Complex& z2) noexcept {
        return !(z1 == z2); 
    }

    constexpr bool operator!=(const Complex& z, const double alpha) noexcept {
        return !(z == alpha); 
    }

    constexpr bool operator!=(const double alpha, const Complex& z) noexcept {
        return !(z == alpha); 
    }

    constexpr std::ostream& operator<< (std::ostream& ostream, const Complex& z) noexcept {
        char sign = std::signbit(z.getY()) ? '-' : '+'; 
        ostream << z.getX() << " " << sign << " " << std::abs(z.getY()) << "j"; 
        return ostream;
    } 

    constexpr std::istream& operator>> (std::istream& istream, Complex& z) noexcept {
        std::string zStr;

        std::string xStr; 
        std::string signStr;
        std::string yStr;

        istream >> xStr >> signStr >> yStr;


        yStr.pop_back();
        if(!isdigit(yStr[yStr.length() - 1])) yStr.pop_back();
        double x = std::stoi(xStr);
        double y = std::stoi(yStr); 
        y = (signStr == "+") ? y : -y;

        z = Complex(x, y);
        
        //z = x + 
        return istream;
    } 

    // add cin op >>. 

    constexpr Complex operator "" _j(const unsigned long long y) noexcept {
        return Complex(0, y);
    }

    constexpr Complex operator "" _j(const long double y) noexcept {
        return Complex(0, y);
    }

    constexpr double real(const Complex& z) noexcept {
        return z.getX();
    }

    constexpr double im(const Complex& z) noexcept {
        return z.getY(); 
    }

    constexpr double mod(const Complex& z) noexcept {
        return std::sqrt(z.getX() * z.getX() + z.getY() * z.getY());
    }

    constexpr double arg(const Complex& z) noexcept {
        return std::atan2(z.getY(), z.getX()); 
    }

    constexpr double norm(const Complex& z) noexcept {
        return z.getX() * z.getX() + z.getY() * z.getY(); 
    }

    constexpr Complex conj(const Complex& z) noexcept {
        return Complex(z.getX(), -z.getY()); 
    }

    constexpr double dot(Complex& z1, Complex& z2) noexcept {
        return  z1.getX() * z2.getX() + z1.getY() * z2.getY(); 
    }

    constexpr Complex proj(const Complex& z) noexcept {
        if(std::isinf(z.getX()) || std::isinf(z.getY())){
            return std::signbit(z.getY()) ? Complex(std::numeric_limits<double>::infinity(), -0.0) : Complex(std::numeric_limits<double>::infinity(), 0);
        }
        return z;
    }

    constexpr Complex polar(const double r, const double theta) noexcept {
        // Consider that polar is one of:
        // r * e^(i * theta)
        // rcos(theta) + i * rsin(theta)
        return Complex(r * std::cos(theta), r * std::sin(theta));
    }

    constexpr Complex exp(const Complex& z) noexcept {
        return std::exp(z.getX()) * Complex(std::cos(z.getY()), std::sin(z.getY()));
    }

    constexpr Complex log(const Complex& z) noexcept {
        return Complex(std::log(mod(z)), arg(z));
    }

    // Will be switched to constexpr by C++26. 
    inline Complex log10(const Complex& z) noexcept {
        return log(z) / std::log(10.0);
    }

    // Will be switched to constexpr by C++26. 
    inline Complex log2(const Complex& z) noexcept {
        return log(z) / std::log(2.0);
    }

    constexpr Complex pow(const Complex& z1, const Complex& z2) noexcept {
        return exp(z2 * log(z1)); 
    }

    constexpr Complex pow(const Complex& z, const double alpha) noexcept {
        return exp(alpha * log(z)); 
    }

    constexpr Complex pow(const double alpha, const Complex& z) noexcept {
        return exp(z * std::log(alpha)); 
    }

    constexpr Complex sqrt(const Complex& z) noexcept {
        double angle = arg(z) / 2; 
        double modZ = mod(z); 

        return std::sqrt(modZ) * Complex(std::cos(angle), std::sin(angle));
    }

    constexpr Complex sin(const Complex& z) noexcept {
        return Complex(std::sin(z.getX()) * std::cosh(z.getY()), std::cos(z.getX()) * std::sinh(z.getY()));
    }

    constexpr Complex cos(const Complex& z) noexcept {
        return Complex(std::cos(z.getX()) * std::cosh(z.getY()), -std::sin(z.getX()) * std::sinh(z.getY()));
    }

    constexpr Complex tan(const Complex& z) noexcept {
        //return Complex(std::tan(z.getX()), std::tanh(z.getY())) / Complex(1, -std::tan(z.getX()) * std::tanh(z.getY())); 
        return sin(z) / cos(z);
    }

    constexpr Complex asin(const Complex& z) noexcept {
        // arcsin(z) = (1/i)ln(iz + sqrt(1 - z^2))
        // arcsin(z) = -i * ln(iz + sqrt(1 - z^2))
        return  -1_j * log(1_j * z + sqrt(1 - z * z)); 
    }

    constexpr Complex acos(const Complex& z) noexcept {
        // return  0.5 * M_PI + 1_j * log(1_j * z + sqrt(1 - z * z)); 
        return  -1_j * log(z + sqrt(z * z - 1));
    }

    constexpr Complex atan(const Complex& z) noexcept {
        return  0.5_j * (log(1 - 1_j * z) - log(1 + 1_j * z)); 
    }

    constexpr Complex sinh(const Complex& z) noexcept {
        return Complex(std::sinh(z.getX()) * std::cos(z.getY()), std::cosh(z.getX()) * std::sin(z.getY())); 
    }

    constexpr Complex cosh(const Complex& z) noexcept {
        return Complex(std::cosh(z.getX()) * std::cos(z.getY()), std::sinh(z.getX()) * std::sin(z.getY())); 
    }

    constexpr Complex tanh(const Complex& z) noexcept {
        // return Complex(std::tanh(z.getX()), std::tan(z.getY())) / Complex(1, std::tanh(z.getX()) * std::tan(z.getY())); 
        return sinh(z) / cosh(z);
    }

    constexpr Complex asinh(const Complex& z) noexcept {
        return  log(z + sqrt(z * z + 1)); 
    }

    constexpr Complex acosh(const Complex& z) noexcept {
        return  log(z + sqrt(z * z - 1)); 
    }

    constexpr Complex atanh(const Complex& z) noexcept {
        return  0.5 * log((1 + z)/(1 - z)); 
    }

    constexpr Complex sgn(const Complex& z) noexcept {
        return z / mod(z); 
    }

    constexpr std::vector<Complex> rootsOfUnity(const int q) noexcept {
        std::vector<Complex> roots(q); 
        Complex omega = exp(1_j * 2 * M_PI / double(q)); 
        // for(int k = 0; k < q; k++){
        //     roots[k] = std::pow(mod(z), 1.0/double(q)) * exp(1_j * arg(z) / double(q)) * pow(omega, k);  
        // }
        // return roots; 

        for(int k = 0; k < q; k++){
            roots[k] = pow(omega, k);  
        }
        return roots; 
    }

    constexpr bool isNan(const Complex& z) noexcept {
        return (std::isnan(z.getX()) || std::isnan(z.getY())) ? true : false; 
    }
}

#endif // COMPLEX_HPP