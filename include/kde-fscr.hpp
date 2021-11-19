#ifndef FSCR_KDE_HPP
#define FSCR_KDE_HPP

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include <type_traits>

#include "kde-kernels.hpp"

namespace fscr
{ 
  class KDE
  { 
    public:
    enum class Bandwith { Scott, Silverman, Custom };

    private:
    inline static double scott_h(const double stdev, const double n) {
      return 1.06 * stdev * std::pow(n, -0.2);
    }
    template<typename T>
    inline static double silverman_h(std::vector<T> data, const double stdev) {
      const size_t n = data.size();
      
      const size_t mid_idx = n / 2;
      size_t Q1_idx, Q3_idx;
      if ((data.size() % 2) == 0) {
        Q1_idx = mid_idx / 2;
        Q3_idx = mid_idx + Q1_idx;
      } else {
        Q1_idx = mid_idx / 2;
        Q3_idx = mid_idx + Q1_idx + 1;
      }
      std::nth_element(data.begin(), data.begin() + Q1_idx, data.end());
      auto Q1 = data[Q1_idx];
      std::nth_element(data.begin(), data.begin() + Q3_idx, data.end());
      auto Q3 = data[Q3_idx];

      double IQR = static_cast<double>(Q3 - Q1);
      
      double IQR_div = IQR / 1.34;
      double A = std::min(stdev, IQR_div);
      return 0.9 * A * std::pow(n, -0.2);
    }

    template<typename T, typename U, typename F>
    static std::vector<double> pdf(const std::vector<T>& data, const std::vector<U>& x_domain, F &&kernel, Bandwith bandwith_type, double bandwith) {
      static_assert(std::is_arithmetic<T>(), "Data types can only be arithmetic (integral or floating-point type");
      static_assert(std::is_arithmetic<U>(), "X domain types can only be arithmetic (integral or floating-point type");
      if (data.size() == 0) {
        std::cerr << "fscr::KDE::pdf() - WARNING: empty data! (1st arg)" << std::endl;
        return std::vector<double>{};
      }
      if (x_domain.size() == 0) {
        std::cerr << "fscr::KDE::pdf() - WARNING: empty x_domain! (2nd arg)" << std::endl;
        return std::vector<double>{};
      }
      const double n = static_cast<double>(data.size());

      double sum = std::accumulate(data.begin(), data.end(), 0.0);
      double mean = sum / n;
      double accum = 0.0;
      for (const auto d: data) {
        accum += (d - mean) * (d - mean);
      }
      double stdev = std::sqrt(accum / (n - 1));

      const auto minmax_itr = std::minmax_element(data.begin(), data.end());
      const double min_val = static_cast<double>(*(minmax_itr.first));
      const double max_val = static_cast<double>(*(minmax_itr.second));
      
      if (bandwith_type == Bandwith::Scott) {
        bandwith = scott_h(stdev, n);
      } else if (bandwith_type == Bandwith::Silverman) {
        bandwith = silverman_h(data, stdev);
      }

      std::vector<double> y_pdf;
      y_pdf.reserve(x_domain.size());

      const double one_nh = 1.0 / (n * bandwith);
      for (const auto x: x_domain) {
        double kernel_sum = std::accumulate(data.begin(), data.end(), 0.0, 
          [bandwith, x, &kernel](double partialSum, T xi) {
            return partialSum + kernel((x - xi) / bandwith);
          }
        );

        y_pdf.push_back(one_nh * kernel_sum);
      }

      return y_pdf;
    }
    
    public:
    /**
     * @brief PDF with kernel parameter - pdf(data, x_domain, kernel)
     */
    template<typename T, typename F>
    static std::vector<double> pdf(const std::vector<T>& data, const std::vector<T>& x_domain, F &&kernel) {
      return pdf(data, x_domain, kernel, Bandwith::Scott, -1.0);
    }

    /**
     * @brief PDF with bandwith selection algorithm parameter - pdf(data, x_domain, bandwith_type)
     */
    template<typename T>
    static std::vector<double> pdf(const std::vector<T>& data, const std::vector<T>& x_domain, Bandwith bandwith_type=Bandwith::Scott) {
      return pdf(data, x_domain, GaussianKernel, bandwith_type, -1.0);
    }
    
    /**
     * @brief PDF with custom bandwith value - pdf(data, x_domain, bandwith_val)
     */
    template<typename T>
    static std::vector<double> pdf(const std::vector<T>& data, const std::vector<T>& x_domain, double bandwith_val) {
      return pdf(data, x_domain, GaussianKernel, Bandwith::Custom, bandwith_val);
    }
    
    /**
     * @brief PDF with kernel parameter and bandwith selection algorithm parameter - pdf(data, x_domain, kernel, bandwith_type)
     */
    template<typename T, typename F>
    static std::vector<double> pdf(const std::vector<T>& data, const std::vector<T>& x_domain, F &&kernel, Bandwith bandwith_type) {
      assert(bandwith_type != Bandwith::Custom);
      return pdf(data, x_domain, kernel, bandwith_type, -1.0);
    }
    
    /**
     * @brief PDF with kernel parameter and custom bandwith - pdf(data, x_domain, kernel, bandwith_val)
     */
    template<typename T, typename F>
    static std::vector<double> pdf(const std::vector<T>& data, const std::vector<T>& x_domain, F &&kernel, double bandwith_val) {
      return pdf(data, x_domain, kernel, Bandwith::Custom, bandwith_val);
    }
    
  };
}

#endif // FSCR_KDE_HPP