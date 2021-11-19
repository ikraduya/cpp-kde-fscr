#ifndef FSCR_KDE_KERNELS_HPP
#define FSCR_KDE_KERNELS_HPP

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32)

#define M_PI   3.14159265358979323846
#define M_PI_2 1.57079632679489661923
#define M_PI_4 0.785398163397448309616
#define M_2_PI 0.636619772367581343076

#else
#define _USE_MATH_DEFINES
#endif
#include <cmath>


namespace fscr
{
  struct {
    inline double operator()(const double x) {
      return 1.0 / std::sqrt(2 * M_PI) * std::exp(-0.5 * x * x);
    }
  } GaussianKernel;

  struct {
    inline double operator()(const double x) {
      return std::abs(x) <= 1.0 ? 0.5 : 0.0;
    }
  } BoxCarKernel; // or Uniform (rectangular window)
  
  struct {
    inline double operator()(const double x) {
      const double abs_x = std::abs(x);
      return abs_x <= 1.0 ? 1.0 - abs_x : 0.0;
    }
  } TriangularKernel;
  
  struct {
    inline double operator()(const double x) {
      return std::abs(x) <= 1.0 ? 0.75 * (1 - x * x) : 0.0;
    }
  } EpanechnikovKernel;
  
  struct {
    inline double operator()(const double x) {
      const double temp = 1 - (x * x);
      return std::abs(x) <= 1.0 ? 0.9375 * (temp * temp) : 0.0;
    }
  } QuarticKernel;
  
  struct {
    inline double operator()(const double x) {
      const double temp = 1 - (x * x);
      return std::abs(x) <= 1.0 ? 1.09375 * (temp * temp * temp) : 0.0;
    }
  } TriweightKernel;
  
  struct {
    inline double operator()(const double x) {
      const double abs_x = std::abs(x);
      const double temp = 1 - (abs_x * abs_x * abs_x);
      return abs_x <= 1.0 ? 0.86419753086 * (temp * temp * temp) : 0.0;
    }
  } TricubeKernel;
  
  struct {
    inline double operator()(const double x) {
      return std::abs(x) <= 1.0 ? M_PI_4 * std::cos(M_PI_2 * x) : 0.0;
    }
  } CosineKernel;
  
  struct {
    inline double operator()(const double x) {
      return 1.0 / (std::exp(x) + 2 + std::exp(-x));
    }
  } LogisticKernel;
  
  struct {
    inline double operator()(const double x) {
      return M_2_PI * (1.0 / (std::exp(x) + std::exp(-x)));
    }
  } SigmoidFunctionKernel;
} 

#endif  // FSCR_KDE_KERNELS_HPP