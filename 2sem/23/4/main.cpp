#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

double a(double x, double y) { return (x * x + y * y + 1.) / 10.; }

double a_x(double x, double y) { return x / 5.; }

double a_y(double x, double y) { return y / 5.; }

double r(double x, double y) {
  return sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5));
}

double f(double x, double y) {
  double r_ = r(x, y);
  double answ = 0.;
  if (r_ < 0.2) {
    answ = exp(-1. / (1. - r_ * r_));
  }
  return answ;
}

std::vector<std::vector<double>> create_computational_map(size_t n) {
  std::vector<std::vector<double>> comp(n + 1);
  for (size_t i = 0; i < n + 1; i++) {
    for (size_t j = 0; j < n + 1; j++) {
      comp.at(i).push_back(0.);
    }
  }
  return comp;
}

void apply_boundary_condition(std::vector<std::vector<double>> &comp,
                              size_t n) {
  for (size_t i = 0; i < n + 1; i++) {
    comp.at(i).at(0) = 0.;
    comp.at(i).at(n) = 0.;
  }
  for (size_t j = 0; j < n + 1; j++) {
    comp.at(0).at(j) = 0.;
    comp.at(n).at(j) = 0.;
  }
}

void method(std::vector<std::vector<double>> &map_0,
            std::vector<std::vector<double>> &map_1, std::vector<double> &error,
            size_t n, double n_t, double h) {
  for (size_t t = 1; t < n_t; t++) {
    std::cout << "Current iteration: " << t << std::endl;

    apply_boundary_condition(map_1, n);

    double error_sum = 0.;

    for (size_t i = 1; i < n; i++) {
      for (size_t j = 1; j < n; j++) {
        map_1.at(i).at(j) =
            h * h * f(i * h, j * h) / (4. * a(i * h, j * h)) +
            0.25 * (map_1.at(i - 1).at(j) + map_0.at(i + 1).at(j) +
                    map_1.at(i).at(j - 1) + map_1.at(i).at(j + 1)) +
            +h / (4. * a(i * h, j * h)) *
                (a_x(i * h, j * h) *
                     (map_0.at(i + 1).at(j) - map_1.at(i - 1).at(j)) +
                 a_y(i * h, j * h) *
                     (map_0.at(i).at(j + 1) - map_1.at(i).at(j - 1)));
        error_sum += abs(map_1.at(i).at(j) - map_0.at(i).at(j));
      }
    }
    error.at(t) = error_sum;
    for (size_t i = 0; i < n + 1; i++) {
      for (size_t j = 0; j < n + 1; j++) {
        map_0.at(i).at(j) = map_1.at(i).at(j);
      }
    }
  }
}

void print_data(std::vector<std::vector<double>> &comp, std::string file_name,
                size_t n, double h) {
  std::ofstream file;
  file.open(file_name);
  for (size_t i = 0; i < n + 1; i++) {
    for (size_t j = 0; j < n + 1; j++) {
      file << i * h << " " << j * h << " " << comp.at(i).at(j) << std::endl;
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

int main(int argc, char *argv[]) {
  const size_t n = 20;
  double n_t = 0;
  double h = 1. / n;

  modf(6. * n * n / (M_PI * M_PI), &n_t);
  n_t *= 3.;

  std::cout << "Number of iterations: " << n_t << std::endl;

  auto comp_0 = create_computational_map(n);
  auto comp_1 = create_computational_map(n);

  std::vector<double> error(n_t);

  apply_boundary_condition(comp_0, n);

  method(comp_0, comp_1, error, n, n_t, h);

  print_data(comp_0, "data.dat", n, h);
  print_error(error, "error.dat", n_t);

  return 0;
}