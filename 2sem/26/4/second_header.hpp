#pragma once


#include "blaze/Math.h"

double f(double x, double y) {
  return -32. * (x * (1. - x) + y * (1. - y));
}

template <size_t n_>
void print_grid(const blaze::StaticMatrix<double, n_, n_, blaze::rowMajor>
                    &computational_grid,
                std::string file_name, double h) {
  std::ofstream file;
  file.open(file_name);
  for (size_t x = 0UL; x < n_; x++) {
    for (size_t y = 0UL; y < n_; y++) {
      file << x * h << " " << y * h << " " << computational_grid.at(x, y)
           << std::endl;
    }
  }
  file.close();
}

void print_error(std::vector<double> &error, std::string file_name,
                 double n_t) {
  std::ofstream error_file;
  error_file.open(file_name);
  for (size_t t = 0; t < n_t; t++) {
    if (t < 2) {
      error_file << t << " " << error.at(t) << " " << 0. << std::endl;
    } else {
      error_file << t << " " << error.at(t) << " "
                 << error.at(t) / (1. - error.at(t - 1) / error.at(t - 2))
                 << std::endl;
    }
  }
  error_file.close();
}