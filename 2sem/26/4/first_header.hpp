#pragma once

#include <cmath>
#include <string>

#include "blaze/Math.h"

double a(int k) { return 0.5 + 0.1 * k; }

double h(double x, double t, int k) {
  return 2. * a(k) / pow(a(k) * (x + t) + 2., 2.);
}

double u_0_t(double t) { return 0.; }

double u_1_t(double t, int k) { return 1. / (a(k) * (t + 1.) + 2.); }

// u_x_0
double f(double x, int k) { return x / (a(k) * x + 2.); }

// u_x_0_t
double g(double x, int k) { return -a(k) * x / pow(a(k) * x + 2., 2.); }

template <size_t nt_, size_t nx_>
void print_grid(const blaze::StaticMatrix<double, nt_, nx_, blaze::rowMajor>
                    &computational_grid,
                std::string file_name, double ht, double hx) {
  std::ofstream file;
  file.open(file_name);
  for (size_t t = 0UL; t < nt_; t++) {
    for (size_t x = 0; x < nx_; x++) {
      file << t * ht << " " << x * hx << " " << computational_grid.at(t, x)
           << std::endl;
    }
  }
  file.close();
}

double analytical_solution(double x, double t, int k) {
  return x / (a(k) * (x + t) + 2.);
}
