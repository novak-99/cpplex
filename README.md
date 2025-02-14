# cpplex

<p align="center">
  <kbd><img src="https://github.com/novak-99/cpplex/blob/main/logo.png?raw=true"/></kbd>
</p>

![C++2b](https://img.shields.io/badge/C++-2023%20(2b)-blue.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Issues](https://img.shields.io/badge/issues-0%20open-red)](https://github.com/novak-99/cpplex/issues)
[![Docs](https://img.shields.io/badge/docs-rtd-green)](https://cpplex.readthedocs.io/en/latest/)

Complex analysis in C++ is slow and lacks a lot of features.

Cpplex makes operations with complex numbers a lot faster and adds much more functionality!

Whether you need tools for scientific computing, signal processing, or numerical analysis, cpplex is the right library for you!

## Usage  

### The Cpplex Complex Datatype

Begin by including the header files of the content you wish to use:

```cpp
#include <Complex.hpp>
#include <Special.hpp>
#include <iostream>
int main() {

  return 0;
}
```

Cpplex allows you to use standard C++ complex literals, as well as our own:

```cpp
#include <Complex.hpp>
#include <Special.hpp>
#include <iostream>

using namespace std::complex_literals; // for C++'s literals

int main() {
  Complex z = 5 + 5_j;
  Complex w = 5 + 5i;

  return 0;
}
```

And you can call functions on the datatype by doing the following:

```cpp
#include <Complex.hpp>
#include <Special.hpp>
#include <iostream>
int main() {
  Complex z = 5 + 5_j;
  Complex gammaZ = gamma(z);

  return 0;
}
```

### Functions in Cpplex

Various functions in cpplex, including the derivative, integral, and continuous entropy functions, require you to use functions of complex inputs. Here is how cpplex handles them.

You can either create complex functions by using lambdas or C++ functions:

```cpp
#include <Complex.hpp>
#include <iostream>

Complex f(Complex z) {
    return z*z; 
}

int main() {

  auto g = [](Complex z) { return z*z };

  return 0;
}
```

Then you can call the relevant function by passing your complex function as a function pointer:

```cpp
#include <Complex.hpp>
#include <NumericalAnalysis.hpp>
#include <iostream>

int main() {
  Complex z = 1 + 1_j;
  auto g = [](Complex z) { return z*z };
  Complex dgdz = derivative(g, z); // Derivative of g evaluated at z.

  return 0;
}
```

### Complex Sequences in Cpplex

Other functions, including discrete transforms and discrete entropy functions, will require complex sequences as input. 

Complex sequences are currently supported in cpplex by using the ```std::vector<Complex>``` datatype, but they will be later swapped out for a vectorized cpplex vector type. 

Complex sequences can be implemented in cpplex by writing the following code:

```cpp
#include <Complex.hpp>
#include <FFT.hpp>
#include <iostream>

int main() {
  std::vector<Complex> X = {1 + 1_j, 2, 3 + 5_j};
  std::vector<Complex> Y = fft(X);

  return 0;
}
```

### Compiling Code

Cpplex is a header only library, meaning no source files or shared object files are required for compilation.

To compile code while using cpplex, remember to include the library's main directory and to compile with C++23:

```g++ main.cpp -I cpplex -o main.o -std=c++2b```

## Documentation

Tutorials for all cpplex modules and detailed documentation for every function are available on the Read the Docs page [here](https://cpplex.readthedocs.io/en/latest/).

## Benchmarks

Cpplex is faster than the standard C++ complex library -- MUCH faster.

For our benchmarks, complex division and exponentiation are **9 times faster**, the complex square root is **10 times faster**, and inverse trig functions are up to **47 times faster**. The full list of benchmarks is available [here](https://cpplex.readthedocs.io/en/latest/Benchmarks/Complex.html).

One of the main reasons for this is that C++'s primitive multiplication and division operators use NaN checking. However, in most complex use cases this is unnecessary and makes programs much slower.

For example, consider the following standard C++ code and the cpplex equivalent for a large sum to approximate the Riemann zeta function at its first zero:
```cpp
void stdZetaSum(){
    const int N = 1e+9;
    std::complex<double> s = 0.5 + 14.1347251417346937904572519835625i;
    std::complex<double> zetaTerm = std::complex<double>(1)/(std::complex<double>(1) - pow(2, std::complex<double>(1) - s));


    std::complex<double> sum = 0; 
    for(int n = 1; n <= N; n++){
        std::complex<double> etaTerm = std::pow(-1, n-1) * 1 / (pow(n, s));

        sum += etaTerm * zetaTerm;
    }
    std::cout << sum << "\n";
}

void cpplexZetaSum(){
    const int N = 1e+9;
    Complex s = 0.5 + 14.1347251417346937904572519835625_j;
    Complex zetaTerm = 1/(1 - pow(2, 1 - s));

    Complex sum = 0; 
    for(int n = 1; n <= N; n++){
        Complex etaTerm = pow(-1, n-1) * 1 / (pow(n, s));

        sum += etaTerm * zetaTerm;
    }
    std::cout << sum << "\n";
}
```
Whereas the former implementation takes an average of 46.788 seconds to run, the cpplex equivalent takes only **18.258 seconds**!

## Future Plans 

### 1. Better approximations 

The Lanczos approximation for the Gamma function and 15-point Gauss-Kronrod quadrature for integration are both good approximation methods, but they both could be better. I hope to swap these out in the future.

### 2. Complex linear algebra module

Currently, cpplex uses the `std::vector<Complex>` type for vectors and does not currently support a complex matrix type. I plan to add an optimized n-d array type in the future.

### 3. Bivariate distributions

Cpplex currently assumes independence for its complex distributions, which can be an erronous assumption. I plan to change the structure of distributions so that they can support correlation between the real and complex components, if it exists.

## Acknowledgements

While the various sources I used for each function can be found in the documentation, in this section, I would like to acknowledge both the [cppreference documentation page](https://en.cppreference.com/w/) and the [scipy documentation page](https://docs.scipy.org/doc/scipy/reference/index.html#scipy-api). Both of these pages helped me see which functions are typically implemented in complex/scientific computing libraries, sources on how to implement them, and naming conventions.
