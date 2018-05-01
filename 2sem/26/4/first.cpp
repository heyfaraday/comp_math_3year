#include <fstream>
#include <iostream>

#include "first_header.hpp"

#include "blaze/Math.h"

int main(int argc, char* argv[]) {
  blaze::setNumThreads(4);

  constexpr int Nx = 200;
  constexpr int Nt = 200;

  constexpr double hx = 1. / Nx;
  constexpr double ht = 1. / Nt;

  constexpr double alpha = ht / hx;

  constexpr int k = 1;

  blaze::StaticMatrix<double, Nt + 1, Nx + 1, blaze::rowMajor> computation_grid;

  // null layer
  for (size_t i = 0UL; i < Nx + 1; ++i) {
    computation_grid.at(0, i) = f(i * hx, k);
  }

  // first layer
  computation_grid.at(1, 0) = u_0_t(ht);
  computation_grid.at(1, Nx) = u_1_t(ht, k);
  for (size_t i = 1UL; i < Nx; ++i) {
    computation_grid.at(1, i) =
        ht * g(i * hx, k) + (1. - alpha * alpha) * computation_grid.at(0, i) +
        0.5 * alpha * alpha *
            (computation_grid.at(0, i - 1) + computation_grid.at(0, i + 1)) +
        0.5 * ht * ht * h(i * hx, 0., k);
  }

  for (size_t t = 1UL; t < Nt; ++t) {
    computation_grid.at(t + 1, 0) = u_0_t(ht * t + ht);
    computation_grid.at(t + 1, Nx) = u_1_t(ht * t + ht, k);
    for (size_t i = 1UL; i < Nx; ++i) {
      computation_grid.at(t + 1, i) =
          -computation_grid.at(t - 1, i) +
          2. * (1. - alpha * alpha) * computation_grid.at(t, i) +
          alpha * alpha *
              (computation_grid.at(t, i - 1) + computation_grid.at(t, i + 1)) +
          ht * ht * h(i * hx, t * ht, k);
    }
  }

  blaze::StaticMatrix<double, Nt + 1, Nx + 1, blaze::rowMajor> analytical_grid;

  for (size_t t = 0; t < Nt + 1; ++t) {
    for (size_t i = 0; i < Nx + 1; ++i) {
      analytical_grid.at(t, i) = analytical_solution(i * hx, t * ht, k);
    }
  }

  blaze::StaticMatrix<double, Nt + 1, Nx + 1, blaze::rowMajor> error_grid;
  for (size_t t = 0; t < Nt + 1; ++t) {
    for (size_t i = 0; i < Nx + 1; ++i) {
      error_grid.at(t, i) = computation_grid.at(t, i) - analytical_grid(t, i);
    }
  }

  // std::cout << computation_grid << std::endl;

  print_grid(computation_grid, "bin/data.dat", ht, hx);
  print_grid(analytical_grid, "bin/analytical_solution.dat", ht, hx);
  print_grid(error_grid, "bin/error.dat", ht, hx);

  return 0;
}