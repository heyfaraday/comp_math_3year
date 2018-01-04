#include <stdio.h>
#include <vector>

#include <cstdlib>
#include <iostream>
#include <math.h>

// You may use this function to access matrix elements, like:
// A(i,j)  ->  array[IND(i,j,N)]  where A is an N-by-N matrix, column-order, arr
// is a N*N vector.
inline int IND(int i, int j, int N) { return i + N * j; }

// You may use this function to read data from ASCII files
void ReadData(const char *fname, std::vector<double> *pV) {
  char buf[1024];
  FILE *fid = fopen(fname, "r");
  if (fid == 0) {
    printf("File not found: %s\n", fname);
    return;
  }

  int line = 0;
  float val;
  while (fgets(buf, 1024, fid) != NULL) {
    sscanf(buf, "%e", &val);
    pV->push_back(val);
    ++line;
  }

  fclose(fid);
}

// You may use this function to write out data to ASCII files
void WriteData(const char *fname, const std::vector<double> &v) {
  FILE *fid = fopen(fname, "w");
  if (fid == 0) {
    printf("Cannot create file: %s\n", fname);
    return;
  }

  for (int n = 0; n < v.size(); ++n) {
    float val = v[n];
    fprintf(fid, "%.6e\n", val);
  }

  fclose(fid);
}

// You may use this function to print data on screen
void PrintData(const std::vector<double> &v) {
  for (int n = 0; n < v.size(); ++n) {
    float val = v[n];
    printf("%.6e\n", val);
  }
}

double f(const std::vector<double> &x) {
  return exp(x.at(0) + 3 * x.at(1) - 0.1) + exp(x.at(0) - 3 * x.at(1) - 0.1) +
         exp(-x.at(0) - 0.1);
}

std::vector<double> grad_f(const std::vector<double> &x) {

  std::vector<double> grad;
  grad.push_back(0.);
  grad.push_back(0.);

  grad.at(0) = exp(x.at(0) + 3. * x.at(1) - 0.1) +
               exp(x.at(0) - 3. * x.at(1) - 0.1) - exp(-x.at(0) - 0.1);
  grad.at(1) = 3. * exp(x.at(0) + 3. * x.at(1) - 0.1) -
               3. * exp(x.at(0) - 3. * x.at(1) - 0.1);

  return grad;
}

std::vector<double> make_step(const std::vector<double> &x) {

  std::vector<double> grad = grad_f(x);

  std::vector<double> step;
  step.push_back(0.);
  step.push_back(0.);

  double H_11 = exp(x.at(0) + 3. * x.at(1) - 0.1) +
                exp(x.at(0) - 3. * x.at(1) - 0.1) + exp(-x.at(0) - 0.1);
  double H_12 = 3. * exp(x.at(0) + 3. * x.at(1) - 0.1) -
                3. * exp(x.at(0) - 3. * x.at(1) - 0.1);
  double H_21 = 3. * exp(x.at(0) + 3. * x.at(1) - 0.1) -
                3. * exp(x.at(0) - 3. * x.at(1) - 0.1);
  double H_22 = 9. * exp(x.at(0) + 3. * x.at(1) - 0.1) +
                9. * exp(x.at(0) - 3. * x.at(1) - 0.1);

  double det = H_11 * H_22 - H_12 * H_21;

  step.at(0) = -(H_22 * grad.at(0) - H_12 * grad.at(1)) / det;
  step.at(1) = -(-H_21 * grad.at(0) + H_11 * grad.at(1)) / det;

  return step;
}

void do_iter(std::vector<double> &x, double learning_rate) {

  std::vector<double> step = make_step(x);

  x.at(0) += learning_rate * step.at(0);
  x.at(1) += learning_rate * step.at(1);
}

void do_method(std::vector<double> &x, double learning_rate,
               int number_of_iterations, std::vector<double> &x_1_history,
               std::vector<double> &x_2_history,
               std::vector<double> &f_x_history) {

  for (int i = 0; i < number_of_iterations; ++i) {
    do_iter(x, learning_rate);
    x_1_history.push_back(x.at(0));
    x_2_history.push_back(x.at(1));
    f_x_history.push_back(f(x));
  }
}

void do_iter_gradient_descent(std::vector<double> &x, double learning_rate) {

  std::vector<double> step = grad_f(x);

  x.at(0) -= learning_rate * step.at(0);
  x.at(1) -= learning_rate * step.at(1);
}

void do_gradient_descent(std::vector<double> &x, double learning_rate,
                         int number_of_iterations,
                         std::vector<double> &x_1_history,
                         std::vector<double> &x_2_history,
                         std::vector<double> &f_x_history) {

  for (int i = 0; i < number_of_iterations; ++i) {
    do_iter_gradient_descent(x, learning_rate);
    x_1_history.push_back(x.at(0));
    x_2_history.push_back(x.at(1));
    f_x_history.push_back(f(x));
  }
}

double Fletcher_Reeves_beta(const std::vector<double> &delta_x_n,
                            const std::vector<double> &delta_x_n_1) {
  return (delta_x_n.at(0) * delta_x_n.at(0) +
          delta_x_n.at(1) * delta_x_n.at(1)) /
         (delta_x_n_1.at(0) * delta_x_n_1.at(0) +
          delta_x_n_1.at(1) * delta_x_n_1.at(1));
}

void do_iter_conj_gradient_method(std::vector<double> &x,
                                  std::vector<double> &previous_s,
                                  std::vector<double> &previous_step) {

  std::vector<double> step = grad_f(x);
  step.at(0) = -step.at(0);
  step.at(1) = -step.at(1);

  double beta = Fletcher_Reeves_beta(step, previous_step);

  std::vector<double> current_s;
  current_s.push_back(0);
  current_s.push_back(0);

  current_s.at(0) = step.at(0) + beta * previous_s.at(0);
  current_s.at(1) = step.at(1) + beta * previous_s.at(1);

  double alpha = 0.003;

  // 1d Newton method
  for (int i = 0; i < 5; ++i) {
    alpha -= (exp(x.at(0) + 3 * x.at(1) - 0.1 +
                  alpha * (current_s.at(0) + 3. * current_s.at(1))) *
                  (current_s.at(0) + 3. * current_s.at(1)) +
              exp(x.at(0) - 3 * x.at(1) - 0.1 +
                  alpha * (current_s.at(0) - 3. * current_s.at(1))) *
                  (current_s.at(0) - 3. * current_s.at(1)) +
              exp(-x.at(0) - 0.1 + alpha * (-current_s.at(0))) *
                  (-current_s.at(0))) /
             (exp(x.at(0) + 3 * x.at(1) - 0.1 +
                  alpha * (current_s.at(0) + 3. * current_s.at(1))) *
                  (current_s.at(0) + 3. * current_s.at(1)) *
                  (current_s.at(0) + 3. * current_s.at(1)) +
              exp(x.at(0) - 3 * x.at(1) - 0.1 +
                  alpha * (current_s.at(0) - 3. * current_s.at(1))) *
                  (current_s.at(0) - 3. * current_s.at(1)) *
                  (current_s.at(0) - 3. * current_s.at(1)) +
              exp(-x.at(0) - 0.1 + alpha * (-current_s.at(0))) *
                  (-current_s.at(0)) * (-current_s.at(0)));
  }

  //  std::cout << alpha << std::endl;

  x.at(0) += alpha * current_s.at(0);
  x.at(1) += alpha * current_s.at(1);

  previous_s = current_s;
  previous_step = step;
}

void do_conj_gradient_method(std::vector<double> &x, int number_of_iterations,
                             std::vector<double> &x_1_history,
                             std::vector<double> &x_2_history,
                             std::vector<double> &f_x_history) {

  std::vector<double> s = grad_f(x);
  s.at(0) = -s.at(0);
  s.at(1) = -s.at(1);

  std::vector<double> step = s;

  do_iter_gradient_descent(x, 0.00001);

  for (int i = 0; i < number_of_iterations; ++i) {
    do_iter_conj_gradient_method(x, s, step);
    x_1_history.push_back(x.at(0));
    x_2_history.push_back(x.at(1));
    f_x_history.push_back(f(x));
  }
}

int main(int argc, char *argv[]) {

  std::vector<double> x;
  x.push_back(-1.);
  x.push_back(2.);

  double learning_rate = 0.01;
  int number_of_iterations = 10;

  std::vector<double> x_1_history;
  std::vector<double> x_2_history;
  std::vector<double> f_x_history;

//    do_method(x, learning_rate, number_of_iterations, x_1_history,
//    x_2_history,
//  			f_x_history);

//    do_gradient_descent(x, learning_rate, number_of_iterations, x_1_history,
//    x_2_history,
//              f_x_history);

  do_conj_gradient_method(x, number_of_iterations, x_1_history, x_2_history,
                          f_x_history);

  WriteData("x_1_2.txt", x_1_history);
  WriteData("x_2_2.txt", x_2_history);
  WriteData("f_x_2.txt", f_x_history);

  return 0;
}