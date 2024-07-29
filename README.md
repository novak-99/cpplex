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

## Benchmarks

## Future Plans 

## Acknowledgements
