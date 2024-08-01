# cpplex

Complex analysis in C++ is slow, lackluster, and boring.

Cpplex makes operations with complex numbers a lot faster, and also adds much more functionality!

Regardless of whether you need tools for scientific computing, signal processing, or numerical analysis, cpplex is the right library for you!

## Usage 

Cpplex is a header only library, meaning no source files or shared object files are required for compilation. 

### The Cpplex Complex Datatype

Begin by including the header files of the content you wish to use:

```cpp
#include <Complex/Complex.hpp>
#include <Special/Special.hpp>
#include <iostream>
int main() {

  return 0;
}
```

Cpplex allows you to use standard C++ complex literals, as well as our own:

```cpp
#include <Complex/Complex.hpp>
#include <Special/Special.hpp>
#include <iostream>

using namespace std::complex_literals; // for C++'s literals

int main() {
  Complex z = 5 + 5_j;
  Complex w = 5 + 5i;

  return 0;
}
```

And you can call functions using this type like so:

```cpp
#include <Complex/Complex.hpp>
#include <Special/Special.hpp>
#include <iostream>
int main() {
  Complex z = 5 + 5_j;
  Complex gammaZ = gamma(z);

  return 0;
}
```

### Functions in Cpplex

Various functions in cpplex, including the derivative, integral, and continuous entropy functions, require you to use functions of complex inputs. Here's how cpplex handles them.

You can either create complex functions by using lambdas or C++ functions:

```cpp
#include <Complex/Complex.hpp>
#include <iostream>

Complex f(Complex z) {
    return z*z; 
}

int main() {

  auto g = [](Complex z) { return z*z };

  return 0;
}
```

Then you can call the relevant function by passing your Complex function as a function pointer:

```cpp
#include <Complex/Complex.hpp>
#include <NumericalAnalysis/NumericalAnalysis.hpp>
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

And so, complex sequences can be implemented in cpplex as such:

```cpp
#include <Complex/Complex.hpp>
#include <FFT/FFT.hpp>
#include <iostream>

int main() {
  std::vector<Complex> X = {1 + 1_j, 2, 3 + 5_j};
  std::vector<Complex> Y = fft(X);

  return 0;
}
```

## Documentation

Tutorials for all cpplex modules and detailed documentation for every function is available on the read the docs page [here](https://cpplex.readthedocs.io/en/latest/).

## Benchmarks

Cpplex is faster than the standard C++ complex library -- MUCH faster.

For our benchmarks, complex division and exponentiation are 9 times faster, the complex square root is 9 times faster, and inverse trig functions are up to 47 times faster. The full list of benchmarks is available [here](https://cpplex.readthedocs.io/en/latest/Benchmarks/Complex.html).

One of the main reasons for this is that C++'s primitive multiplication and division operators use NaN checking. In most complex use cases however, this is unnecessary, and makes programs much slower.

## Future Plans 

### 1. Better approximations 

The Lanczos approximation for the Gamma function and 15-point Gauss-Kronrod quadrature for integration are both good approximation methods, but they both could be better. I hope to swap these out in the future.

### 2. Complex linear algebra module

Currently, cpplex uses the `std::vector<Complex>` type for vectors and does not currently support a complex matrix type. I plan to add an optimized n-d array type in the future.

### 3. Bivariate Distributions

Cpplex currently assumes independence for its complex distributions, which can be an erronous assumption. I plan to change the structure of distributions so that they can support correlation between the real and complex components, if it exists.

## Acknowledgements

While the various sources I used for each function can be found in the documentation, here, I would like to specially thank both the [cppreference documentation page](cppreference.com) and the [scipy documentation page](https://docs.scipy.org/doc/scipy/reference/index.html#scipy-api). Both of these pages helped me see which functions are typically implemented in complex/scientific computing libraries, sources on how to implement them, and naming conventions.
