#include "Complex.hpp"
#include <cmath>
#include <limits>
#ifdef __x86_64__
 #include <immintrin.h>
#else   
 #include <arm_neon.h>
#endif


namespace cpplex {

    // Complex::Complex(){
    //     x = 0; 
    //     y = 0; 
    // }

    // Complex::Complex(const double x, const double y){
    //     this->x = x; 
    //     this->y = y; 
    // }

    // void Complex::setX(const double x){
    //     this->x = x; 
    // }

    // void Complex::setY(const double y){
    //     this->y = y; 
    // }

    // double Complex::getX() const{
    //     return x; 
    // }

    // double Complex::getY() const{
    //     return y; 
    // }

    // double Complex::real(){
    //     return x; 
    // }

    // double Complex::im(){
    //     return y; 
    // }

    // Complex& Complex::operator+=(const Complex& z){
    //     x += z.getX(); 
    //     y += z.getY(); 
    //     return *this; 
    // }

    // Complex& Complex::operator-=(const Complex& z){
    //     x -= z.getX(); 
    //     y -= z.getY(); 
    //     return *this; 
    // }

    // Complex& Complex::operator*=(const Complex& z){

    //     // (a * c - b * d, a * d + c * b)

    //     // double a = z1.getX(); 
    //     // double b = z1.getY();
    //     // double c = z2.getX();  
    //     // double d = z2.getY(); 

    //     x = x * z.getX() - y * z.getY(); 
    //     y = x * z.getY() + y * z.getX(); 

    //     return *this; 
    // }

    // Complex& Complex::operator/=(const Complex& z){
    //     // double a = z1.getX(); 
    //     // double b = z1.getY();
    //     // double c = z2.getX();  
    //     // double d = z2.getY(); 

    //     // double den = c * c + d * d; 
    //     // return Complex((a * c + b * d)/den,  (b * c - a * d)/den); 

    //     double den = z.getX() * z.getX() + z.getY() * z.getY(); 
    //     x = (x * z.getX() + y * z.getY()) / den; 
    //     y = (y * z.getX() - x * z.getY()) / den; 
    //     return *this; 
    // }

    // Complex& Complex::operator+=(const double alpha){
    //     x += alpha; 
    //     return *this;  
    // }

    // Complex& Complex::operator-=(const double alpha){
    //     x -= alpha; 
    //     return *this; 
    // }

    // Complex& Complex::operator*=(const double alpha){
    //     x *= alpha; 
    //     y *= alpha;
    //     return *this; 
    // }

    // Complex& Complex::operator/=(const double alpha){
    //     x /= alpha; 
    //     y /= alpha;
    //     return *this; 
    // }

    // double real(const Complex& z){
    //     return z.getX();
    // }

    // double im(const Complex& z){
    //     return z.getY(); 
    // }

    // double mod(const Complex& z){
    //     return std::sqrt(z.getX() * z.getX() + z.getY() * z.getY());
    // }

    // double arg(const Complex& z){
    //     return std::atan2(z.getY(), z.getX()); 
    // }

    // double norm(const Complex& z){
    //     return mod(z) * mod(z);
    // }

    // Complex conj(const Complex& z){
    //     return Complex(z.getX(), -z.getY()); 
    // }

    // double dot(Complex& z1, Complex& z2){
    //     return  z1.getX() * z2.getX() + z1.getY() * z2.getY(); 
    // }

    // Complex proj(const Complex& z){
    //     if(std::isinf(z.getX()) || std::isinf(z.getY())){
    //         return std::signbit(z.getY()) ? Complex(std::numeric_limits<double>::infinity(), -0) : Complex(std::numeric_limits<double>::infinity(), 0);
    //     }
    //     return z;
    // }

    // Complex polar(const double r, const double theta){
    //     // Consider that polar is one of:
    //     // r * e^(i * theta)
    //     // rcos(theta) + i * rsin(theta)
    //     return Complex(r * std::cos(theta), r * std::sin(theta));
    // }

    // Complex operator+(const Complex& z){
    //     return z; 
    // }

    // Complex operator-(const Complex& z){
    //     return Complex(-z.getX(), -z.getY());
    // }

    // Complex operator+(const Complex& z1, const Complex& z2){
    //     return Complex(z1.getX() + z2.getX(), z1.getY() + z2.getY()); 
    // }

    // Complex operator-(const Complex& z1, const Complex& z2){
    //     return Complex(z1.getX() - z2.getX(), z1.getY() - z2.getY()); 
    // }

    // Complex operator*(const Complex& z1, const Complex& z2){
    //     double a = z1.getX(); 
    //     double b = z1.getY();
    //     double c = z2.getX();  
    //     double d = z2.getY(); 

    //     return Complex(a * c - b * d, a * d + c * b); 
    // }

    // Complex operator/(const Complex& z1, const Complex& z2){
    //     double a = z1.getX(); 
    //     double b = z1.getY();
    //     double c = z2.getX();  
    //     double d = z2.getY(); 

    //     double den = c * c + d * d; 
    //     return Complex((a * c + b * d)/den,  (b * c - a * d)/den); 
    // }

    // Complex operator+(const Complex& z, const double alpha){
    //     return Complex(z.getX() + alpha, z.getY()); 
    // }

    // Complex operator-(const Complex& z, const double alpha){
    //     return Complex(z.getX() - alpha, z.getY()); 
    // }

    // Complex operator*(const Complex& z, const double alpha){
    //     return Complex(z.getX() * alpha, z.getY() * alpha); 
    // }

    // Complex operator/(const Complex& z, const double alpha){
    //     return Complex(z.getX() / alpha, z.getY() / alpha); 
    // }

    // Complex operator+(const double alpha, const Complex& z){
    //     return Complex(z.getX() + alpha, z.getY()); 
    // }

    // Complex operator-(const double alpha, const Complex& z){
    //     return Complex(alpha - z.getX(), -z.getY()); 
    // }

    // Complex operator*(const double alpha, const Complex& z){
    //     return Complex(z.getX() * alpha, z.getY() * alpha); 
    // }

    // Complex operator/(const double alpha, const Complex& z){
    //     double den = norm(z); 
    //     return Complex(alpha * z.getX() / den, -alpha * z.getY() / den); 
    // }

    // bool operator==(const Complex& z1, const Complex& z2){
    //     return true ? (z1.getX() == z2.getX() && z1.getY() == z2.getY()) : false; 
    // }

    // bool operator==(const Complex& z, const double alpha){
    //     return true ? (z.getX() == alpha && z.getY() == 0) : false; 
    // }

    // bool operator==(const double alpha, const Complex& z){
    //     return true ? (z.getX() == alpha && z.getY() == 0) : false; 
    // }

    // bool operator!=(const Complex& z1, const Complex& z2){
    //     return !(z1 == z2); 
    // }

    // bool operator!=(const Complex& z, const double alpha){
    //     return !(z == alpha); 
    // }

    // bool operator!=(const double alpha, const Complex& z){
    //     return !(z == alpha); 
    // }

    // Complex operator "" _j(const unsigned long long y){
    //     return Complex(0, y);
    // }

    // Complex operator "" _j(const long double y){
    //     return Complex(0, y);
    // }

    // std::ostream& operator<< (std::ostream& ostream, const Complex& z){
    //     char sign = std::signbit(z.getY()) ? '-' : '+'; 
    //     ostream << z.getX() << " " << sign << " " << std::abs(z.getY()) << "j"; 
    //     return ostream;
    // }


}