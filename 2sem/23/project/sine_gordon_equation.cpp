#include <cmath>
#include <fstream>
#include <vector>

double analytical_solution(double x, double t);

double f(double x);

double g(double x);

double kink(double x, double t) {
  double x0 = -5.;
  double c = 0.3;
  return 4. * atan(exp((x - x0 - c * t) / sqrt(1. - c * c)));
}

// double antikink(double x, double t) {
//   double x0 = -5.;
//   double c = 0.3;
//   return 4. * atan(exp(-(x - x0 - c * t) / sqrt(1. - c * c)));
// }

// double kink_kink_collision(double x, double t) {
//   double c = 0.6;
//   return 4. *
//          atan(c * sinh(x / sqrt(1. - c * c)) / cosh(c * t / sqrt(1. - c *
//          c)));
// }

// double kink_antikink_collision(double x, double t) {
//   double c = 0.6;
//   return 4. *
//          atan(sinh(c * t / sqrt(1. - c * c)) / cosh(x / sqrt(1. - c * c)) /
//          c);
// }

// double breather(double x, double t) {
//   double omega = 0.4;
//   return 4. * atan(sqrt(1. - omega * omega) * sin(omega * t) / omega /
//                    cosh(sqrt(1. - omega * omega) * x));
// }

double kink_derivative_at_0(double x) {
  double x0 = -5.;
  double c = 0.3;
  return -2. * (c / sqrt(1. - c * c)) *
         (1. / cosh((x - x0) / sqrt(1. - c * c)));
}

// double antikink_derivative_at_0(double x) {
//   double x0 = -5.;
//   double c = 0.3;
//   return 2. * (c / sqrt(1. - c * c)) * (1. / cosh((x - x0) / sqrt(1. - c * c)));
// }

// double breather_derivative_at_0(double x) {
//   double omega = 0.4;
//   return 4. * sqrt(1. - omega * omega) / cosh(x * sqrt(1. - omega * omega));
// }

// double kink_kink_collision_derivative_at_0(double x) {
//   return 0.;
// }

// double kink_antikink_collision_derivative_at_0(double x) {
//   double c = 0.6;
//   return 4. / sqrt(1. - c * c) / cosh(x / sqrt(1. - c * c));
// }

double analytical_solution(double x, double t) {
  return kink(x, t);
  // return antikink(x, t);
  // return kink_kink_collision(x, t);
  // return kink_antikink_collision(x, t);
  // return breather(x, t);
}

double f(double x) {
  return kink(x, 0);
  // return antikink(x, 0);
  // return kink_kink_collision(x, 0);
  // return kink_antikink_collision(x, 0);
  // return breather(x, 0);
  return 0.;
}

double g(double x) {
  return kink_derivative_at_0(x);
  // return antikink_derivative_at_0(x);
  // return kink_kink_collision_derivative_at_0(x);
  // return kink_antikink_collision_derivative_at_0(x);
  // return breather_derivative_at_0(x);
}

int main(int argc, char *argv[]) {
  size_t nx = 400;
  double lower_bound_x = -20.;
  double upper_bound_x = 20.;
  double hx = (upper_bound_x - lower_bound_x) / nx;
  size_t nt = 1800;
  double upper_bound_t = 90;
  double ht = (upper_bound_t - 0.) / nt;

  double alpha = ht / hx;

  std::vector<std::vector<double>> comp(nt + 1);
  for (size_t t = 0; t < nt + 1; t++) {
    for (size_t i = 0; i < nx + 1; i++) {
      comp.at(t).push_back(0.);
    }
  }

  // null layer
  for (size_t i = 0; i < nx + 1; i++) {
    comp.at(0).at(i) = f(i * hx + lower_bound_x);
  }

  // first layer
  for (size_t i = 1; i < nx; i++) {
    comp.at(1).at(i) =
        ht * g(i * hx + lower_bound_x) +
        (1. - alpha * alpha) * comp.at(0).at(i) +
        0.5 * alpha * alpha * (comp.at(0).at(i - 1) + comp.at(0).at(i + 1)) -
        0.5 * ht * ht * sin(comp.at(0).at(i));
  }
  // comp.at(1).at(0) =
  //     ht * g(lower_bound_x) + (1. - alpha * alpha) * comp.at(0).at(0) +
  //     0.5 * alpha * alpha * (comp.at(0).at(1) + comp.at(0).at(1)) -
  //     0.5 * ht * ht * sin(comp.at(0).at(0));
  // comp.at(1).at(nx) =
  //     ht * g(upper_bound_x) +
  //     (1. - alpha * alpha) * comp.at(0).at(nx) +
  //     0.5 * alpha * alpha * (comp.at(0).at(nx - 1) + comp.at(0).at(nx - 1)) -
  //     0.5 * ht * ht * sin(comp.at(0).at(nx));

  comp.at(1).at(0) = analytical_solution(lower_bound_x, ht);
  comp.at(1).at(nx) = analytical_solution(upper_bound_x, ht);

  // for (size_t i = 0; i < nx + 1; i++) {
  //   comp.at(1).at(i) = analytical_solution(i * hx + lower_bound_x, ht);
  // }

  // next layers
  for (size_t t = 1; t < nt; t++) {
    for (size_t i = 1; i < nx; i++) {
      comp.at(t + 1).at(i) =
          -comp.at(t - 1).at(i) + 2. * (1. - alpha * alpha) * comp.at(t).at(i) +
          alpha * alpha * (comp.at(t).at(i - 1) + comp.at(t).at(i + 1)) -
          ht * ht * sin(comp.at(t).at(i));
    }
    // comp.at(t + 1).at(0) =
    //     -comp.at(t - 1).at(0) + 2. * (1. - alpha * alpha) * comp.at(t).at(0)
    //     + alpha * alpha * (comp.at(t).at(1) + comp.at(t).at(1)) - ht * ht *
    //     sin(comp.at(t).at(0));
    // comp.at(t + 1).at(nx) =
    //     -comp.at(t - 1).at(nx) + 2. * (1. - alpha * alpha) *
    //     comp.at(t).at(nx) + alpha * alpha * (comp.at(t).at(nx - 1) +
    //     comp.at(t).at(nx - 1)) - ht * ht * sin(comp.at(t).at(nx));
    comp.at(t + 1).at(0) = analytical_solution(lower_bound_x, ht * (t + 1));
    comp.at(t + 1).at(nx) = analytical_solution(upper_bound_x, ht * (t + 1));
  }

  std::ofstream analytical_solution_file;
  analytical_solution_file.open("analytical_solution.dat");
  for (size_t t = 0; t < nt + 1; t++) {
    for (size_t x = 0; x < nx + 1; x++) {
      analytical_solution_file
          << t * ht << " " << x * hx + lower_bound_x << " "
          << analytical_solution(x * hx + lower_bound_x, t * ht) << std::endl;
    }
  }
  analytical_solution_file.close();

  std::ofstream numerical_solution_file;
  numerical_solution_file.open("numerical_solution.dat");
  for (size_t t = 0; t < nt + 1; t++) {
    for (size_t x = 0; x < nx + 1; x++) {
      numerical_solution_file << t * ht << " " << x * hx + lower_bound_x << " "
                              << comp.at(t).at(x) << std::endl;
    }
  }
  numerical_solution_file.close();

  std::ofstream error_file;
  error_file.open("error.dat");
  for (size_t t = 0; t < nt + 1; t++) {
    for (size_t x = 0; x < nx + 1; x++) {
      error_file << t * ht << " " << x * hx + lower_bound_x << " "
                 << comp.at(t).at(x) -
                        analytical_solution(x * hx + lower_bound_x, t * ht)
                 << std::endl;
    }
  }
  error_file.close();

  return 0;
}