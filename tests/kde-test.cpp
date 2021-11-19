#include <gtest/gtest.h>

#include <cmath>

#include "kde-fscr.hpp"
#include "kde-kernels.hpp"

namespace {
  // absolute different for EXPECT_NEAR() floating-point values
  constexpr auto absoluteError = 0.00001;

  template<typename T=double>
  std::vector<T> linspace(double min_val, double max_val, size_t num) {
    std::vector<T> ret(num);
    T interval = (max_val - min_val) / (num - 1);
    for (int i=0; i<num; ++i) {
      ret[i] = min_val + (interval * i);
    }
    return ret;
  }
} //anonymous namespace

TEST(KernelFunc, tc1GaussianKernel) {
  EXPECT_NEAR(fscr::GaussianKernel(0), 0.3989422804, absoluteError);
  EXPECT_NEAR(fscr::GaussianKernel(1), 0.24197072451, absoluteError);
  EXPECT_NEAR(fscr::GaussianKernel(10), 7.6945986e-23, absoluteError);
  EXPECT_NEAR(fscr::GaussianKernel(-1), 0.24197072451, absoluteError);
  EXPECT_NEAR(fscr::GaussianKernel(-10), 7.6945986e-23, absoluteError);
}

TEST(KernelFunc, tc2BoxCarKernel) {
  EXPECT_NEAR(fscr::BoxCarKernel(0), 0.5, absoluteError);
  EXPECT_NEAR(fscr::BoxCarKernel(1), 0.5, absoluteError);
  EXPECT_NEAR(fscr::BoxCarKernel(10), 0, absoluteError);
  EXPECT_NEAR(fscr::BoxCarKernel(-1), 0.5, absoluteError);
  EXPECT_NEAR(fscr::BoxCarKernel(-10), 0, absoluteError);
}

TEST(KernelFunc, tc3TriangularKernel) {
  EXPECT_NEAR(fscr::TriangularKernel(0), 1, absoluteError);
  EXPECT_NEAR(fscr::TriangularKernel(0.5), 0.5, absoluteError);
  EXPECT_NEAR(fscr::TriangularKernel(10), 0, absoluteError);
  EXPECT_NEAR(fscr::TriangularKernel(-0.5), 0.5, absoluteError);
  EXPECT_NEAR(fscr::TriangularKernel(-10), 0, absoluteError);
}

TEST(KDE_PDF, tc1EmptyData) {
  const std::vector<double> series{};
  const size_t num = 10;
  const std::vector<double> x_domain = linspace(-7, 11, num);
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain);

  EXPECT_EQ(y_pdf.size(), 0);
}

TEST(KDE_PDF, tc2EmptyXDomain) {
  const std::vector<double> series{6.2, 5.1, 1.9, -0.4, -1.3, -2.1};
  const std::vector<double> x_domain{};
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain);

  EXPECT_EQ(y_pdf.size(), 0);
}

TEST(KDE_PDF, tc3EmptyDataFloat) {
  const std::vector<float> series{};
  const size_t num = 10;
  const std::vector<float> x_domain = linspace<float>(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain);

  EXPECT_EQ(y_pdf.size(), 0);
}

TEST(KDE_PDF, tc4EmptyXDomainFloat) {
  const std::vector<float> series{6.2f, 5.1f, 1.9f, -0.4f, -1.3f, -2.1f};
  const std::vector<float> x_domain{};
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain);

  EXPECT_EQ(y_pdf.size(), 0);
}

TEST(KDE_PDF_KGaussianBScott, tc1NormalInput) {
  const std::vector<double> series{6.2, 5.1, 1.9, -0.4, -1.3, -2.1};
  const size_t num = 10;
  const std::vector<double> x_domain = linspace(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain);
  const std::vector<double> expectedValue{
    0.007307, 0.028649, 0.065150, 0.090492, 0.087123, 
    0.074772, 0.066339, 0.048622, 0.022996, 0.006330
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KGaussianBScott, tc2NormalInputFloat) {
  const std::vector<float> series{6.2f, 5.1f, 1.9f, -0.4f, -1.3f, -2.1f};
  const size_t num = 10;
  const std::vector<float> x_domain = linspace<float>(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain);
  const std::vector<double> expectedValue{
    0.007307, 0.028649, 0.065150, 0.090492, 0.087123, 
    0.074772, 0.066339, 0.048622, 0.022996, 0.006330
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}


TEST(KDE_PDF_KGaussianBSilverman, tc1NormalInput) {
  const std::vector<double> series{6.2, 5.1, 1.9, -0.4, -1.3, -2.1};
  const size_t num = 10;
  const std::vector<double> x_domain = linspace(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::KDE::Bandwith::Silverman);
  const std::vector<double> expectedValue{
    0.003684, 0.023164, 0.068024, 0.100008, 0.088391, 
    0.071691, 0.069916, 0.051552, 0.019584, 0.003428
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KGaussianBSilverman, tc2NormalInputFloat) {
  const std::vector<float> series{6.2f, 5.1f, 1.9f, -0.4f, -1.3f, -2.1f};
  const size_t num = 10;
  const std::vector<float> x_domain = linspace<float>(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::KDE::Bandwith::Silverman);
  const std::vector<double> expectedValue{
    0.003684, 0.023164, 0.068024, 0.100008, 0.088391, 
    0.071691, 0.069916, 0.051552, 0.019584, 0.003428
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KGaussianBCustom, tc1NormalInput) {
  const std::vector<double> series{6.2, 5.1, 1.9, -0.4, -1.3, -2.1};
  const size_t num = 10;
  const std::vector<double> x_domain = linspace(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, std::sqrt(2.25));
  const std::vector<double> expectedValue{
    0.000249, 0.009358, 0.070429, 0.125095, 0.085787,
    0.059328, 0.081730, 0.058461, 0.009273, 0.000284
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KGaussianBCustom, tc2NormalInputFloat) {
  const std::vector<float> series{6.2f, 5.1f, 1.9f, -0.4f, -1.3f, -2.1f};
  const size_t num = 10;
  const std::vector<float> x_domain = linspace<float>(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, std::sqrt(2.25));
  const std::vector<double> expectedValue{
    0.000249, 0.009358, 0.070429, 0.125095, 0.085787,
    0.059328, 0.081730, 0.058461, 0.009273, 0.000284
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KBoxCarBScott, tc1NormalInput) {
  const std::vector<double> series{6.2, 5.1, 1.9, -0.4, -1.3, -2.1};
  const size_t num = 10;
  const std::vector<double> x_domain = linspace(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::BoxCarKernel);
  const std::vector<double> expectedValue{
    0.000000, 0.000000, 0.065172, 0.097758, 0.097758,
    0.065172, 0.065172, 0.065172, 0.000000, 0.000000
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KBoxCarBScott, tc2NormalInputFloat) {
  const std::vector<float> series{6.2f, 5.1f, 1.9f, -0.4f, -1.3f, -2.1f};
  const size_t num = 10;
  const std::vector<float> x_domain = linspace<float>(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::BoxCarKernel);
  const std::vector<double> expectedValue{
    0.000000, 0.000000, 0.065172, 0.097758, 0.097758,
    0.065172, 0.065172, 0.065172, 0.000000, 0.000000
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KBoxCarBSilverman, tc1NormalInput) {
  const std::vector<double> series{6.2, 5.1, 1.9, -0.4, -1.3, -2.1};
  const size_t num = 10;
  const std::vector<double> x_domain = linspace(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::BoxCarKernel, fscr::KDE::Bandwith::Silverman);
  const std::vector<double> expectedValue{
    0.000000, 0.000000, 0.076758, 0.115137, 0.076758,
    0.076758, 0.076758, 0.076758, 0.000000, 0.000000
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KBoxCarBSilverman, tc2NormalInputFloat) {
  const std::vector<float> series{6.2f, 5.1f, 1.9f, -0.4f, -1.3f, -2.1f};
  const size_t num = 10;
  const std::vector<float> x_domain = linspace<float>(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::BoxCarKernel, fscr::KDE::Bandwith::Silverman);
  const std::vector<double> expectedValue{
    0.000000, 0.000000, 0.076758, 0.115137, 0.076758,
    0.076758, 0.076758, 0.076758, 0.000000, 0.000000
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KBoxCarBCustom, tc1NormalInput) {
  const std::vector<double> series{6.2, 5.1, 1.9, -0.4, -1.3, -2.1};
  const size_t num = 10;
  const std::vector<double> x_domain = linspace(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::BoxCarKernel, std::sqrt(2.25));

  const std::vector<double> expectedValue{
    0.000000, 0.000000, 0.055556, 0.166667, 0.111111,
    0.055556, 0.111111, 0.055556, 0.000000, 0.000000
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KBoxCarBCustom, tc2NormalInputFloat) {
  const std::vector<float> series{6.2f, 5.1f, 1.9f, -0.4f, -1.3f, -2.1f};
  const size_t num = 10;
  const std::vector<float> x_domain = linspace<float>(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::BoxCarKernel, std::sqrt(2.25));
  const std::vector<double> expectedValue{
    0.000000, 0.000000, 0.055556, 0.166667, 0.111111,
    0.055556, 0.111111, 0.055556, 0.000000, 0.000000
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KTriangularBScott, tc1NormalInput) {
  const std::vector<double> series{6.2, 5.1, 1.9, -0.4, -1.3, -2.1};
  const size_t num = 10;
  const std::vector<double> x_domain = linspace(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::TriangularKernel);

  const std::vector<double> expectedValue{
    0.000000, 0.000000, 0.064085, 0.144547, 0.078288,
    0.048794, 0.097214, 0.061536, 0.000000, 0.000000
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KTriangularBScott, tc2NormalInputFloat) {
  const std::vector<float> series{6.2f, 5.1f, 1.9f, -0.4f, -1.3f, -2.1f};
  const size_t num = 10;
  const std::vector<float> x_domain = linspace<float>(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::TriangularKernel);
  const std::vector<double> expectedValue{
    0.000000, 0.000000, 0.064085, 0.144547, 0.078288,
    0.048794, 0.097214, 0.061536, 0.000000, 0.000000
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KTriangularBSilverman, tc1NormalInput) {
  const std::vector<double> series{6.2, 5.1, 1.9, -0.4, -1.3, -2.1};
  const size_t num = 10;
  const std::vector<double> x_domain = linspace(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::TriangularKernel, fscr::KDE::Bandwith::Silverman);

  const std::vector<double> expectedValue{
    0.000000, 0.000000, 0.061604, 0.159572, 0.072209,
    0.040394, 0.107560, 0.058069, 0.000000, 0.000000
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KTriangularBSilverman, tc2NormalInputFloat) {
  const std::vector<float> series{6.2f, 5.1f, 1.9f, -0.4f, -1.3f, -2.1f};
  const size_t num = 10;
  const std::vector<float> x_domain = linspace<float>(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::TriangularKernel, fscr::KDE::Bandwith::Silverman);
  const std::vector<double> expectedValue{
    0.000000, 0.000000, 0.061604, 0.159572, 0.072209,
    0.040394, 0.107560, 0.058069, 0.000000, 0.000000
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KTriangularBCustom, tc1NormalInput) {
  const std::vector<double> series{6.2, 5.1, 1.9, -0.4, -1.3, -2.1};
  const size_t num = 10;
  const std::vector<double> x_domain = linspace(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::TriangularKernel, std::sqrt(2.25));

  const std::vector<double> expectedValue{
    0.000000, 0.000000, 0.044444, 0.185185, 0.051852,
    0.029630, 0.125926, 0.051852, 0.000000, 0.000000
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}

TEST(KDE_PDF_KTriangularBCustom, tc2NormalInputFloat) {
  const std::vector<float> series{6.2f, 5.1f, 1.9f, -0.4f, -1.3f, -2.1f};
  const size_t num = 10;
  const std::vector<float> x_domain = linspace<float>(-7, 11, num);
  
  const std::vector<double> y_pdf = fscr::KDE::pdf(series, x_domain, fscr::TriangularKernel, std::sqrt(2.25));
  const std::vector<double> expectedValue{
    0.000000, 0.000000, 0.044444, 0.185185, 0.051852,
    0.029630, 0.125926, 0.051852, 0.000000, 0.000000
  };

  EXPECT_EQ(y_pdf.size(), expectedValue.size());
  for (auto i=0; i<y_pdf.size(); ++i) {
    EXPECT_NEAR(y_pdf[i], expectedValue[i], absoluteError);
  }
}
