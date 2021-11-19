#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "cpp-kde-fscr/kde-fscr.hpp"

std::vector<double> linspace(double min_val, double max_val, size_t num) {
  std::vector<double> ret(num);
  double interval = (max_val - min_val) / (num - 1);
  for (int i=0; i<num; ++i) {
    ret[i] = min_val + (interval * i);
  }
  return ret;
}

void print_vec(const std::vector<double> &vec) {
  for (const auto v: vec) {
    std::cout << v << ", ";
  }
  std::cout << '\n' << std::endl;
}

int main(int argc, char const *argv[]) {
  std::cout << std::setprecision(6);

  std::vector<double> series = {6.2, 5.1, 1.9, -0.4, -1.3, -2.1};
  const int nums = 10;
  std::vector<double> x_domain = linspace(-7, 11, nums);

  std::cout << "x_domain: " << std::endl;
  print_vec(x_domain);

  std::vector<double> y_pdf_kgauss_hscott = fscr::KDE::pdf(series, x_domain);
  std::cout << "PDF with Gaussian kernel and Scott bandwith selection: " << std::endl;
  print_vec(y_pdf_kgauss_hscott);
  
  std::vector<double> y_pdf_ktriangular_hcustom = fscr::KDE::pdf(series, x_domain, fscr::GaussianKernel, std::sqrt(2.25));
  std::cout << "PDF with Triangular kernel and user specified bandwith: " << std::endl;
  print_vec(y_pdf_ktriangular_hcustom);
  
  // Epanechikov (parabolic)
  auto customKernel = [](const double x) -> double {
    return std::abs(x) <= 1.0 ? (0.75 * (1 - (x * x))) : 0.0;
  };
  std::vector<double> y_pdf_kcustom_hsilverman = fscr::KDE::pdf(series, x_domain, customKernel, fscr::KDE::Bandwith::Silverman);
  std::cout << "PDF with Custom kernel and Silverman bandwith selection: " << std::endl;
  print_vec(y_pdf_kcustom_hsilverman);

  return 0;
}