#include <cmath>
#include <fstream>
#include <iostream>

#include "second_header.hpp"

#include "blaze/Math.h"

int main(int argc, char *argv[]) {
  blaze::setNumThreads(4);

  constexpr size_t N = 100;
  constexpr double h = 1. / N;
  double l = M_PI * M_PI * 2.;
  double L = 8. / pow(h, 2.) - l;
  double rho = 1. - 2. * l / L;
  double rho_mpn = 1. - 2. * sqrt(l / L);
  double Nt = -6. / log(rho);
  double Nt_mpn = -6. / log(rho_mpn);

  modf(Nt * 4., &Nt);
  modf(Nt_mpn * 4., &Nt_mpn);

  std::cout << Nt << std::endl;
  std::cout << Nt_mpn << std::endl;

  // double tau_opt_mpi = 2. / (L + l);
  double tau_opt_mpn = 1. / sqrt(l * L);

  blaze::StaticMatrix<double, N + 1, N + 1, blaze::rowMajor> computation_grid_0;
  blaze::StaticMatrix<double, N + 1, N + 1, blaze::rowMajor> computation_grid_1;

  std::vector<double> error(Nt_mpn + 1);

  // null layer
  for (size_t i = 0UL; i < N + 1; ++i) {
    computation_grid_0.at(i, 0) = 0.;
    computation_grid_0.at(i, N) = 0.;
  }
  for (size_t j = 0UL; j < N + 1; ++j) {
    computation_grid_0.at(0, j) = 0.;
    computation_grid_0.at(N, j) = 0.;
  }

  // for (size_t t = 1; t < Nt + 1; ++t) {
  //   // std::cout << "Current iteration: " << t << std::endl;

  //   double error_sum = 0;

  //   for (size_t i = 1UL; i < N; ++i) {
  //     for (size_t j = 1UL; j < N; ++j) {
  //       computation_grid_1.at(i, j) =
  //           computation_grid_0.at(i, j) +
  //           tau_opt_mpi * ((computation_grid_0.at(i - 1, j) +
  //                           computation_grid_0.at(i + 1, j) -
  //                           2. * computation_grid_0.at(i, j)) /
  //                              (h * h) +
  //                          (computation_grid_0.at(i, j - 1) +
  //                           computation_grid_0.at(i, j + 1) -
  //                           2. * computation_grid_0.at(i, j)) /
  //                              (h * h)) -
  //           f(i * h, j * h) * tau_opt_mpi;
  //       error_sum +=
  //           abs(computation_grid_1.at(i, j) - computation_grid_0.at(i, j));
  //     }
  //   }
  //   error.at(t) = error_sum;
  //   for (size_t i = 0UL; i < N + 1; ++i) {
  //     computation_grid_1.at(i, 0) = 0.;
  //     computation_grid_1.at(i, N) = 0.;
  //   }
  //   for (size_t j = 0UL; j < N + 1; ++j) {
  //     computation_grid_1.at(0, j) = 0.;
  //     computation_grid_1.at(N, j) = 0.;
  //   }
  //   computation_grid_0 = computation_grid_1;
  // }

  for (size_t t = 1; t < Nt_mpn + 1; ++t) {
    // std::cout << "Current iteration: " << t << std::endl;

    double error_sum = 0;

    blaze::StaticMatrix<double, N + 1, N + 1, blaze::rowMajor>
        computation_grid_1_2;

    for (size_t j = 1UL; j < N; ++j) {
      computation_grid_1_2.at(0, j) = 0.;
      computation_grid_1_2.at(N, j) = 0.;

      std::vector<double> pn(N + 1);
      std::vector<double> qn(N + 1);

      pn.at(1) = 0.;
      qn.at(1) = 0.;

      for (size_t i = 1; i < N; i++) {
        double an = -1. / (h * h);
        double cn = an;
        double bn = -1. / tau_opt_mpn - 2. / (h * h);
        double dn =
            computation_grid_0.at(i, j) / tau_opt_mpn +
            (computation_grid_0.at(i, j - 1) + computation_grid_0.at(i, j + 1) -
             2. * computation_grid_0.at(i, j)) /
                (h * h) -
            f(i * h, j * h);
        pn.at(i + 1) = cn / (bn - an * pn.at(i));
        qn.at(i + 1) = (an * qn.at(i) - dn) / (bn - an * pn.at(i));
      }
      for (size_t i = N - 1; i > 0; i--) {
        computation_grid_1_2.at(i, j) =
            pn.at(i + 1) * computation_grid_1_2.at(i + 1, j) + qn.at(i + 1);
      }
    }

    for (size_t i = 1UL; i < N; ++i) {
      computation_grid_1.at(i, 0) = 0.;
      computation_grid_1.at(i, N) = 0.;

      std::vector<double> pn(N + 1);
      std::vector<double> qn(N + 1);

      pn.at(1) = 0.;
      qn.at(1) = 0.;

      for (size_t j = 1; j < N; j++) {
        double an = -1. / (h * h);
        double cn = an;
        double bn = -1. / tau_opt_mpn - 2. / (h * h);
        double dn = computation_grid_1_2.at(i, j) / tau_opt_mpn +
                    (computation_grid_1_2.at(i - 1, j) +
                     computation_grid_1_2.at(i + 1, j) -
                     2. * computation_grid_1_2.at(i, j)) /
                        (h * h) -
                    f(i * h, j * h);
        pn.at(j + 1) = cn / (bn - an * pn.at(j));
        qn.at(j + 1) = (an * qn.at(j) - dn) / (bn - an * pn.at(j));
      }
      for (size_t j = N - 1; j > 0; j--) {
        computation_grid_1.at(i, j) =
            pn.at(j + 1) * computation_grid_1.at(i, j + 1) + qn.at(j + 1);
      }
    }

    for (size_t i = 1UL; i < N; ++i) {
      for (size_t j = 1UL; j < N; ++j) {
        error_sum +=
            abs(computation_grid_1.at(i, j) - computation_grid_0.at(i, j));
      }
    }
    error.at(t) = error_sum;
    for (size_t i = 0UL; i < N + 1; ++i) {
      computation_grid_1.at(i, 0) = 0.;
      computation_grid_1.at(i, N) = 0.;
    }
    for (size_t j = 0UL; j < N + 1; ++j) {
      computation_grid_1.at(0, j) = 0.;
      computation_grid_1.at(N, j) = 0.;
    }
    computation_grid_0 = computation_grid_1;
  }


print_grid(computation_grid_0, "bin/data.dat", h);
print_error(error, "bin/error.dat", Nt_mpn + 1);

return 0;
}