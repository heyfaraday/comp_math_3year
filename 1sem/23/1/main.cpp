#include <stdio.h>
#include <vector>

#include <cstdlib>
#include <iostream>
#include <math.h>

using namespace std;

// You may use this function to access matrix elements, like:
// A(i,j)  ->  array[IND(i,j,N)]  where A is an N-by-N matrix, column-order, arr
// is a N*N vector.
inline int IND(int i, int j, int N) { return i + N * j; }

// You may use this function to read data from ASCII files
void ReadData(const char *fname, vector<double> *pV) {
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
void WriteData(const char *fname, const vector<double> &v) {
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
void PrintData(const vector<double> &v) {
  for (int n = 0; n < v.size(); ++n) {
    float val = v[n];
    printf("%.6e\n", val);
  }
}

inline double A(int i, int j, int N, double alpha, double beta) {
  if (i == j) {
    if (i == 1) {
      return 2 + beta;
    } else {
      return 2 + alpha;
    }
  } else if (abs(i - j) == 1) {
    return -1;
  } else {
    return 0;
  }
}

inline double B(int i, int j, int N, double alpha, double beta) {
  if (i == j) {
    return 1;
  } else if (abs(i - j) == 1) {
    if (i == 1) {
      return -1.0 / (2.0 + beta);
    } else {
      return -1.0 / (2.0 + alpha);
    }
  } else {
    return 0;
  }
}

inline double I(int i, int j, int N) {
  if (i == j) {
    return 1;
  } else {
    return 0;
  }
}

vector<double> inv_D_f(vector<double> const &f, double alpha, double beta) {
  int N = static_cast<int>(f.size());
  vector<double> f_new(N);
  for (int i = 0; i < N; i++) {
    f_new.at(i) = f.at(i) / A(i, i, N, alpha, beta);
  }
  return f_new;
}

void make_iter_1(vector<double> &x_next, vector<double> const &x_current,
                 vector<double> const &f, double parameter, double alpha,
                 double beta) {

  int N = static_cast<int>(f.size());
  for (int i = 0; i < N; i++) {

    double sum = 0;
    x_next.at(i) = 0;

    if (i != 0 && i != (N - 1)) {
      for (int j = i - 1; j < i + 2; j++) {
        sum +=
            (I(i + 1, j + 1, N) - parameter * A(i + 1, j + 1, N, alpha, beta)) *
            x_current.at(j);
      }
    } else if (i == 0) {
      for (int j = 0; j < 2; j++) {
        sum +=
            (I(i + 1, j + 1, N) - parameter * A(i + 1, j + 1, N, alpha, beta)) *
            x_current.at(j);
      }
    } else if (i == (N - 1)) {
      for (int j = i - 1; j < N; j++) {
        sum +=
            (I(i + 1, j + 1, N) - parameter * A(i + 1, j + 1, N, alpha, beta)) *
            x_current.at(j);
      }
    }

    x_next.at(i) = sum + parameter * f.at(i);
  }
}

void make_iter_2(vector<double> &x_next, vector<double> const &x_current,
                 vector<double> const &f, double parameter, double alpha,
                 double beta) {

  int N = static_cast<int>(f.size());
  for (int i = 0; i < N; i++) {

    double sum = 0;
    x_next.at(i) = 0;

    if (i != 0 && i != (N - 1)) {
      for (int j = i - 1; j < i + 2; j++) {
        sum +=
            (I(i + 1, j + 1, N) - parameter * B(i + 1, j + 1, N, alpha, beta)) *
            x_current.at(j);
      }
    } else if (i == 0) {
      for (int j = 0; j < 2; j++) {
        sum +=
            (I(i + 1, j + 1, N) - parameter * B(i + 1, j + 1, N, alpha, beta)) *
            x_current.at(j);
      }
    } else if (i == (N - 1)) {
      for (int j = i - 1; j < N; j++) {
        sum +=
            (I(i + 1, j + 1, N) - parameter * B(i + 1, j + 1, N, alpha, beta)) *
            x_current.at(j);
      }
    }
    x_next.at(i) = sum + parameter * f.at(i);
  }
}

double Residual(vector<double> const &x_k, vector<double> const &x_0,
                vector<double> const &f, double alpha, double beta) {

  unsigned long N = x_0.size();
  double r_k = 0.0;
  double r_0 = 0.0;

  for (int i = 0; i < N; i++) {

    double sum_k = 0;
    double sum_0 = 0;

    if (i != 0 && i != (N - 1)) {
      for (int j = i - 1; j < i + 2; j++) {
        sum_k += A(i + 1, j + 1, N, alpha, beta) * x_k.at(j);
        sum_0 += A(i + 1, j + 1, N, alpha, beta) * x_0.at(j);
      }
    } else if (i == 0) {
      for (int j = 0; j < 2; j++) {
        sum_k += A(i + 1, j + 1, N, alpha, beta) * x_k.at(j);
        sum_0 += A(i + 1, j + 1, N, alpha, beta) * x_0.at(j);
      }
    } else if (i == (N - 1)) {
      for (int j = i - 1; j < N; j++) {
        sum_k += A(i + 1, j + 1, N, alpha, beta) * x_k.at(j);
        sum_0 += A(i + 1, j + 1, N, alpha, beta) * x_0.at(j);
      }
    }

    r_k += (f.at(i) - sum_k) * (f.at(i) - sum_k);
    r_0 += (f.at(i) - sum_0) * (f.at(i) - sum_0);
  }

  return sqrt(r_k / r_0);
}

vector<double> Solve_1(vector<double> const &x_0, vector<double> const &f,
                       vector<double> &res, double parameter,
                       int number_of_iterations, double alpha, double beta,
                       double border) {
  unsigned long N = x_0.size();
  vector<double> x_k(N);
  vector<double> x_k_1(N);

  x_k = x_0;
  for (int i_iter = 0; i_iter < number_of_iterations; i_iter++) {
    double res_current = Residual(x_k, x_0, f, alpha, beta);
    if (res_current > border) {
      res.push_back(res_current);
      make_iter_1(x_k_1, x_k, f, parameter, alpha, beta);
      x_k = x_k_1;
    } else {
      break;
    }
  }

  return x_k;
}

vector<double> Solve_2(vector<double> const &x_0, vector<double> const &f,
                       vector<double> &res, double parameter,
                       int number_of_iterations, double alpha, double beta,
                       double border) {
  unsigned long N = x_0.size();
  vector<double> x_k(N);
  vector<double> x_k_1(N);

  x_k = x_0;
  for (int i_iter = 0; i_iter < number_of_iterations; i_iter++) {
    double res_current = Residual(x_k, x_0, f, alpha, beta);
    if (res_current > border) {
      res.push_back(res_current);
      make_iter_2(x_k_1, x_k, inv_D_f(f, alpha, beta), parameter, alpha, beta);
      x_k = x_k_1;
    } else {
      break;
    }
  }

  return x_k;
}

int main(int argc, char *argv[]) {

  printf("Start.\n");

  double alpha = 0.01;
  double beta = 10.0;
  double border = 0.001;

  const int number_of_iterations = 4000;

  vector<double> r;
  ReadData("../r10.txt", &r);

  // Initial x
  vector<double> x_0(r.size());
  for (int i = 0; i < r.size(); i++) {
    x_0.at(i) = 0.;
  }

  vector<double> res_1;
  vector<double> res_2;

  vector<double> x_out_1 =
      Solve_1(x_0, r, res_1, 0.15, number_of_iterations, alpha, beta, border);
  vector<double> x_out_2 =
      Solve_2(x_0, r, res_2, 1.0, number_of_iterations, alpha, beta, border);

  WriteData("x_out_1.txt", x_out_1);
  WriteData("x_out_2.txt", x_out_2);

  WriteData("res_1.txt", res_1);
  WriteData("res_2.txt", res_2);

  return 0;
}
