# cpp-kde-fscr | C++ Kernel Density Estimation Library - From Scratch 

**cpp-kde-fscr** is a header-only C++ library of [Kernel Density Estimation (KDE)](https://en.wikipedia.org/wiki/Kernel_density_estimation) functions created without any dependencies. This library is created From Scratch (fscr), literally. 

In addition to provide a C++ library without dependencies, the **fscr** libraries are my attemp to make an easy-to-learn libraries for people who want to learn, especially, Kernel Density Estimation (KDE) code implementation from scratch.

Features:
- A header-only library of KDE functions with easily changeable kernel
- Written in C++11 format, and is C++11/14/17 compatible
- Released under a permissive, MIT license


### Kernels
Available kernels [(wikipedia)](https://en.wikipedia.org/wiki/Kernel_(statistics)):
- Gaussian / Normal (default)
- Boxcar / Uniform
- Triangular
- Epachenikov / Parabolic
- Quartic / Biweight
- Triweight
- Tricube
- Cosine
- Logistic
- Sigmoid function

### Bandwith Selection Algorithms
- Scott
- Silvermann
- Custom (user input)

## Installation
Simply add the repository as the `git submodule` and/or add the header files to your project using:
``` C++
#include "cpp-kde-fscr/kde-fscr.hpp"
```

> Don't forget to add the installation directory into your include path

### CMake
Or, You can install the library from source using CMake.
``` bash
# clone cpp-kde-fscr from GitHub
git clone https://github.com/ikraduya/cpp-kde-fscr

# create build folder and specify installation location
cd cpp-kde-fscr
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/cpp-kde-fscr/installation/location

# install
cmake --install build

# Build tests (optional)
cmake -S . -B build_tests -DBUILD_TESTS=1
cmake --build build_tests --target tests
cd build_tests && ctest
```

For example, `/cpp-kde-fscr/installation/location` could be `/usr/local`.

If `DCMAKE_INSTALL_PREFIX` is not provided, the default installation location for linux would be `/usr/local`.

### Example Usage with CMake
You can see at `example_project/` directory on how to use it in your project. 

## Author
Ikraduya Edian

## License
MIT
