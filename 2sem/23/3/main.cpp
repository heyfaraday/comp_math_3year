#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

double f(double x, double t) {
  return exp(-3. * t) * sin(M_PI * x) * (-3. + M_PI * M_PI * (1. + pow(x, 4.)));
}

double analytical_solution(double x, double t) {
  return exp(-3. * t) * sin(M_PI * x);
}

// For tridiagonal matrix algorithm
double a_n(double x, double hx) { return -(1. + pow(x, 4.)) / (2. * hx * hx); }

double c_n(double x, double hx) { return -(1. + pow(x, 4.)) / (2. * hx * hx); }

double b_n(double x, double hx, double ht) {
  return -1. / ht - (1. + pow(x, 4.)) / (hx * hx);
}

std::vector<std::vector<double>> create_computational_map(int nx, int nt) {
  // First index for t
  // Second index for x

  std::vector<std::vector<double>> comp(nt + 1);
  for (size_t t = 0; t < nt + 1; t++) {
    for (size_t i = 0; i < nx + 1; i++) {
      comp.at(t).push_back(0.);
    }
  }
  return comp;
}

void apply_initial_condition(std::vector<std::vector<double>> &comp, int nx,
                             double hx) {
  for (size_t i = 0; i < nx + 1; i++) {
    comp.at(0).at(i) = sin(M_PI * i * hx);
  }
}

void print_computational_map(const std::vector<std::vector<double>> &comp,
                             std::string file_name, int nx, int nt, double hx,
                             double ht) {
  std::ofstream file;
  file.open(file_name);
  for (size_t t = 0; t < nt + 1; t++) {
    for (size_t i = 0; i < nx + 1; i++) {
      file << ht * t << " " << i * hx << " " << comp.at(t).at(i) << std::endl;
    }
  }
  file.close();
}

void print_analytical_solution(std::string file_name, int nx, int nt, double hx,
                               double ht) {
  std::ofstream file;
  file.open(file_name);
  for (size_t t = 0; t < nt + 1; t++) {
    for (size_t i = 0; i < nx + 1; i++) {
      file << ht * t << " " << i * hx << " "
           << analytical_solution(i * hx, t * ht) << std::endl;
    }
  }
  file.close();
}

void print_solutions_at_the_middle(
    const std::vector<std::vector<double>> &comp,
    const std::vector<std::vector<double>> &comp_2, std::string file_name,
    int nx, int nt, double hx, double ht) {
  std::vector<double> sol(nt + 1);
  for (size_t t = 0; t < nt + 1; t++) {
    sol.at(t) = analytical_solution(0.5, t * ht);
  }
  std::ofstream file;
  file.open(file_name);
  for (size_t t = 0; t < nt + 1; t++) {
    file << ht * t << " " << sol.at(t) << " " << comp.at(t).at(nx / 2) << " "
         << comp_2.at(t).at(nx / 2) << std::endl;
  }
  file.close();
}

void explicit_scheme(std::vector<std::vector<double>> &comp, int nx, int nt,
                     double hx, double ht) {
  for (size_t t = 0; t < nt; t++) {
    comp.at(t + 1).at(0) = 0.;
    comp.at(t + 1).at(nx) = 0.;

    for (size_t i = 1; i < nx; i++) {
      comp.at(t + 1).at(i) = comp.at(t).at(i) + ht * f(i * hx, t * ht) +
                             ht * (1. + pow(i * hx, 4.)) *
                                 (comp.at(t).at(i - 1) + comp.at(t).at(i + 1) -
                                  2. * comp.at(t).at(i)) /
                                 (hx * hx);
    }
  }
}

void crank_nicolson_scheme(std::vector<std::vector<double>> &comp, int nx,
                           int nt, double hx, double ht) {
  for (size_t t = 0; t < nt; t++) {
    comp.at(t + 1).at(0) = 0.;
    comp.at(t + 1).at(nx) = 0.;

    std::vector<double> pn(nx + 1);
    std::vector<double> qn(nx + 1);

    pn.at(1) = 0.;
    qn.at(1) = 0.;

    // Tridiagonal matrix algorithm
    for (size_t i = 1; i < nx; i++) {
      double dn = f(i * hx, t * ht + ht / 2.) + comp.at(t).at(i) / ht +
                  (1. + pow(i * hx, 4)) * (1. / (2. * hx * hx)) *
                      (comp.at(t).at(i + 1) + comp.at(t).at(i - 1) -
                       2. * comp.at(t).at(i));

      pn.at(i + 1) =
          c_n(i * hx, hx) / (b_n(i * hx, hx, ht) - a_n(i * hx, hx) * pn.at(i));
      qn.at(i + 1) = (a_n(i * hx, hx) * qn.at(i) - dn) /
                     (b_n(i * hx, hx, ht) - a_n(i * hx, hx) * pn.at(i));
    }

    for (size_t i = nx - 1; i > 0; i--) {
      comp.at(t + 1).at(i) =
          pn.at(i + 1) * comp.at(t + 1).at(i + 1) + qn.at(i + 1);
    }
  }
}

int main(int argc, char *argv[]) {
  const int nx = 20;
  const int nt = 2000;

  const double hx = 1. / nx;
  const double ht = 1. / nt;

  auto comp = create_computational_map(nx, nt);
  auto comp_2 = create_computational_map(nx, nt);
  apply_initial_condition(comp, nx, hx);
  apply_initial_condition(comp_2, nx, hx);

  explicit_scheme(comp, nx, nt, hx, ht);
  crank_nicolson_scheme(comp_2, nx, nt, hx, ht);

  print_computational_map(comp, "computational_map_explicit.dat", nx, nt, hx,
                          ht);
  print_computational_map(comp_2, "computational_map_cn.dat", nx, nt, hx, ht);
  print_analytical_solution("analytical_solution.dat", nx, nt, hx, ht);

  print_solutions_at_the_middle(comp, comp_2, "solutions_at_the_middle.dat", nx,
                                nt, hx, ht);

  return 0;
}