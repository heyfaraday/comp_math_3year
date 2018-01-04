#include <stdio.h>
#include <vector>

#include <math.h>
#include <iostream>
#include <cstdlib>

using namespace std;

//You may use this function to access matrix elements, like:
//A(i,j)  ->  array[IND(i,j,N)]  where A is an N-by-N matrix, column-order, arr is a N*N vector.
inline int IND(int i,int j,int N){return i+N*j;}

//You may use this function to read data from ASCII files
void ReadData(const char * fname, vector<double> * pV)
{
    char buf[1024];
    FILE * fid = fopen(fname,"r");
    if(fid ==0)
    {
        printf("File not found: %s\n",fname);
        return;
    }

    int line=0;
    float val;
    while (fgets(buf, 1024, fid) != NULL)
    {
        sscanf(buf,"%e",&val);
        pV->push_back(val);
        ++line;
    }

    fclose(fid);
}


//You may use this function to write out data to ASCII files
void WriteData(const char * fname,const vector<double> & v)
{
    FILE * fid = fopen(fname,"w");
    if(fid ==0)
    {
        printf("Cannot create file: %s\n",fname);
        return;
    }

    for(int n = 0; n < v.size();++n)
    {
        float val = v[n];
        fprintf(fid,"%.6e\n",val);
    }

    fclose(fid);
}

//You may use this function to print data on screen
void PrintData(const vector<double> & v)
{
    for(int n = 0; n < v.size();++n)
    {
        float val = v[n];
        printf("%.6e\n",val);
    }

}

inline double A(int i, int j, int N, double alpha, double beta) {
    if (i == j) {
        if (i == 1) {
            return 2 + beta;
        } else {
            return 2 + alpha;
        }
    } else if (abs(i-j) == 1){
        return -1;
    } else {
        return 0;
    }
}

inline double B(int i, int j, int N, double alpha, double beta) {
    if (i == j) {
        return 1;
    } else if (abs(i-j) == 1){
        if (i == 1) {
            return - 1.0 / (2.0 + beta);
        } else {
            return - 1.0 / (2.0 + alpha);
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

void Mult(vector<double> const& A, vector<double> const& B, vector<double> &answ, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {

            for (int k = 0; k < size; k++) {

                answ.at(IND(i, j, size)) += A.at(IND(i, k, size)) * B.at(IND(k, j, size));
            }
        }
    }
}

void Mult_2(vector<double> const& A, vector<double> const& x, vector<double> &answ, int size) {
    for (int i = 0; i < size; i++) {
        for (int k = 0; k < size; k++) {
            answ.at(i) += A.at(IND(i, k, size)) * x.at(k);
        }
    }
}

void make_iter_1(vector<double> &x_next, vector<double> const &x_current, vector<double> const &f, double parameter, double alpha, double beta) {

    int N = static_cast<int>(f.size());
    for (int i = 0; i < N; i++) {

        double sum = 0;
        x_next.at(i) = 0;

        if (i != 0 && i != (N - 1)) {
            for (int j = i - 1; j < i + 2; j++) {
                sum += (I(i + 1, j + 1, N) - parameter * A(i + 1, j + 1, N, alpha, beta)) * x_current.at(j);
            }
        } else if (i == 0) {
            for (int j = 0; j < 2; j++) {
                sum += (I(i + 1, j + 1, N) - parameter * A(i + 1, j + 1, N, alpha, beta)) * x_current.at(j);
            }
        } else if (i == (N - 1)) {
            for (int j = i - 1; j < N; j++) {
                sum += (I(i + 1, j + 1, N) - parameter * A(i + 1, j + 1, N, alpha, beta)) * x_current.at(j);
            }
        }

        x_next.at(i) = sum + parameter * f.at(i);
    }
}

void make_iter_2(vector<double> &x_next, vector<double> const &x_current, vector<double> const &f, double parameter, double alpha, double beta) {

    int N = static_cast<int>(f.size());
    for (int i = 0; i < N; i++) {

        double sum = 0;
        x_next.at(i) = 0;

        if (i != 0 && i != (N - 1)) {
            for (int j = i - 1; j < i + 2; j++) {
                sum += (I(i + 1, j + 1, N) - parameter * B(i + 1, j + 1, N, alpha, beta)) * x_current.at(j);
            }
        } else if (i == 0) {
            for (int j = 0; j < 2; j++) {
                sum += (I(i + 1, j + 1, N) - parameter * B(i + 1, j + 1, N, alpha, beta)) * x_current.at(j);
            }
        } else if (i == (N - 1)) {
            for (int j = i - 1; j < N; j++) {
                sum += (I(i + 1, j + 1, N) - parameter * B(i + 1, j + 1, N, alpha, beta)) * x_current.at(j);
            }
        }
        x_next.at(i) = sum + parameter * f.at(i);
    }
}

void make_iter_3(vector<double> inv_L_D__U, vector<double> inv_L_D__f, vector<double> &x_next, vector<double> const &x_current, int size) {

    int N = static_cast<int>(x_next.size());

    vector<double> first_term(size);
    Mult_2(inv_L_D__U, x_current, first_term, size);
    vector<double> second_term(size);
    second_term = inv_L_D__f;

    for (int i = 0; i < N; i++) {
        x_next.at(i) = - first_term.at(i) + second_term.at(i);
    }
}

double Residual(vector<double> const &x_k, vector<double> const &x_0, vector<double> const &f, double alpha, double beta) {

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


double Residual_1(vector<double> const& A, vector<double> const &x_k, vector<double> const &x_0, vector<double> const &f) {

    unsigned long N = x_0.size();
    double r_k = 0.0;
    double r_0 = 0.0;

    for (int i = 0; i < N; i++) {

        double sum_k = 0;
        double sum_0 = 0;

        for (int j = 0; j < N; j++) {
            sum_k += A.at(IND(i, j, N)) * x_k.at(j);
            sum_0 += A.at(IND(i, j, N)) * x_0.at(j);
        }

        r_k += (f.at(i) - sum_k) * (f.at(i) - sum_k);
        r_0 += (f.at(i) - sum_0) * (f.at(i) - sum_0);
    }

    return sqrt(r_k / r_0);
}

vector<double> Solve_1(vector<double> const &x_0, vector<double> const& f, double parameter, int number_of_iterations, double alpha, double beta)
{
    unsigned long N = x_0.size();
    vector<double> x_k(N);
    vector<double> x_k_1(N);

    x_k = x_0;
    for (int i_iter = 0; i_iter < number_of_iterations; i_iter++) {
        std::cout << Residual(x_k, x_0, f, alpha, beta) << std::endl;
        make_iter_1(x_k_1, x_k, f, parameter, alpha, beta);
        x_k = x_k_1;
    }

    return x_k;
}

vector<double> Solve_2(vector<double> const &x_0, vector<double> const& f, double parameter, int number_of_iterations, double alpha, double beta)
{
    unsigned long N = x_0.size();
    vector<double> x_k(N);
    vector<double> x_k_1(N);

    x_k = x_0;
    for (int i_iter = 0; i_iter < number_of_iterations; i_iter++) {
        std::cout << Residual(x_k, x_0, f, alpha, beta) << std::endl;
        make_iter_2(x_k_1, x_k, inv_D_f(f, alpha, beta), parameter, alpha, beta);
        x_k = x_k_1;
    }

    return x_k;
}

vector<double> Solve_3(vector<double> const &A,
                        vector<double> const &inv_L_D__U,
                       vector<double> const &inv_L_D__f,
                       vector<double> const &x_0,
                       vector<double> const& f,
                       int number_of_iterations,
                       int size) {

    unsigned long N = x_0.size();
    vector<double> x_k(N);
    vector<double> x_k_1(N);

    x_k = x_0;
    for (int i_iter = 0; i_iter < number_of_iterations; i_iter++) {
//        std::cout << Residual_1(A, x_k, x_0, f) << std::endl;
        make_iter_3(inv_L_D__U, inv_L_D__f, x_k_1, x_k, size);
        x_k = x_k_1;
    }

    return x_k;
}

void Inverse(vector<double> const &A, vector<double> &inv_A, int size) {

    int i = 0, j = 0, k = 0, n = size;
    float **mat = NULL;
    float d = 0.0;

    mat = new float*[2*n];
    for (i = 0; i < 2*n; ++i)
    {
        mat[i] = new float[2*n]();
    }

    for(i = 0; i < n; ++i)
    {
        for(j = 0; j < n; ++j)
        {
            mat[i][j] = A.at(IND(i, j, n));
        }
    }

//    cout << endl << "Input matrix:" << endl;
//    for (i = 0; i < n; ++i)
//    {
//        for (j = 0; j < n; ++j)
//        {
//            cout << mat[i][j] << "\t";
//        }
//        cout << endl;
//    }
//    cout << endl;

//    Initializing Right-hand side to identity matrix
    for(i = 0; i < n; ++i)
    {
        for(j = 0; j < 2*n; ++j)
        {
            if(j == (i+n))
            {
                mat[i][j] = 1;
            }
        }
    }

//    Partial pivoting
    for(i = n; i > 1; --i)
    {
        if(mat[i-1][1] < mat[i][1])
        {
            for(j = 0; j < 2*n; ++j)
            {
                d = mat[i][j];
                mat[i][j] = mat[i-1][j];
                mat[i-1][j] = d;
            }
        }
    }
    cout << endl;

//    Pivoted output
//    cout << "Pivoted output: " << endl;
//    for(i = 0; i < n; ++i)
//    {
//        for(j = 0; j < 2*n; ++j)
//        {
//            cout << mat[i][j] << "\t";
//        }
//        cout << endl;
//    }
//    cout << endl;

//    Reducing To Diagonal Matrix
    for(i = 0; i < n; ++i)
    {
        for(j = 0; j < 2*n; ++j)
        {
            if(j != i)
            {
                d = mat[j][i] / mat[i][i];
                for(k = 0; k < n*2; ++k)
                {
                    mat[j][k] -= mat[i][k]*d;
                }
            }
        }
    }

//    Reducing To Unit Matrix
    for(i = 0; i < n; ++i)
    {
        d = mat[i][i];
        for(j = 0; j < 2*n; ++j)
        {
            mat[i][j] = mat[i][j]/d;
        }
    }

//    Print inverse of the input matrix
//    cout<< "Inverse matrix:" << endl;
//    for(i = 0; i < n; ++i)
//    {
//        for(j = n; j < 2*n; ++j)
//        {
//            cout << mat[i][j] << "\t";
//        }
//        cout << endl;
//    }

    for(i = 0; i < n; ++i)
    {
        for(j = n; j < 2*n; ++j)
        {
           inv_A.at(IND(i, j - n, n)) =  mat[i][j];
        }
    }

    // Deleting the memory allocated
    for (i = 0; i < n; ++i)
    {
        delete[] mat[i];
    }
    delete[] mat;
}

double first_norm(vector<double> const &A, int size) {

    double sum = 0.;
    for (int j = 0; j < size; j++) {
        sum += A[IND(0, j, size)];
    }
    double sum_max = sum;
    for (int i = 1; i < size; i++) {
        sum = 0.;
        for (int j = 0; j < size; j++) {
            sum += A[IND(i, j, size)];
        }
        if (sum > sum_max) {
            sum_max = sum;
        }
    }
    return sum_max;
}

double scalp(vector<double> const &v1, vector<double> const &v2, int size) {
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum += v1[i] * v2[i];
    }
}

void compute_max_and_min_eig(vector<double> const &A, int size, int iter, double &max, double &min) {

    vector<double> u_0(size);
    vector<double> u_k(size);
    vector<double> u_k_1(size);

    for (int i = 0; i < size; i++) {
        u_0.at(i) = 1.;
    }

    u_k = u_0;

    for (int i_iter = 0; i_iter < iter; i_iter++) {
        Mult_2(A, u_k, u_k_1, size);
        u_k = u_k_1;
    }
    max = scalp(u_k_1, u_k, size) / scalp(u_k, u_k, size);

    vector<double> B(size * size);

    for (int i = 0; i < size; i++) {
        B[IND(i, i, size)] = A[IND(i, i, size)] - max;
    }

    u_k = u_0;

    for (int i_iter = 0; i_iter < iter; i_iter++) {
        Mult_2(B, u_k, u_k_1, size);
        u_k = u_k_1;
    }

    min = max - scalp(u_k_1, u_k, size) / scalp(u_k, u_k, size);
}




int main(int argc, char *argv[])
{
    printf("Start.\n");

    const int number_of_iterations = 10;
    const int SIZE = 10;

    vector<double> r(SIZE);
    for (int i = 0; i < SIZE; i++){
        r.at(i) = i + 1;
    }
//    WriteData("../out_r.txt", r);

    vector<double> A(SIZE * SIZE);
    for (int i = 0; i < SIZE; i++){
        A.at(IND(i, i, SIZE)) = 100;
        for (int j = 0; j < i; j++){
            A.at(IND(i, j, SIZE)) = (i + 1.0) / (j + 1.0);
        }
        if (i != SIZE - 1) {
            A.at(IND(i, i + 1, SIZE)) = (i + 1.0) / (i + 2.0);
        }
    }
//    WriteData("../out_A.txt", A);

    vector<double> L_D(SIZE * SIZE);
    for (int i = 0; i < SIZE; i++){
        L_D.at(IND(i, i, SIZE)) = 100;
        for (int j = 0; j < i; j++){
            L_D.at(IND(i, j, SIZE)) = (i + 1.0) / (j + 1.0);
        }
    }

    vector<double> U(SIZE * SIZE);
    for (int i = 0; i < SIZE; i++){
        if (i != SIZE - 1) {
            U.at(IND(i, i + 1, SIZE)) = (i + 1.0) / (i + 2.0);
        }
    }

    vector<double> x_0(r.size());
    for (int i = 0; i < r.size(); i++) {
        x_0.at(i) = 1.;
    }

    vector<double> inv_A(SIZE * SIZE);
    Inverse(A, inv_A, SIZE);
//    WriteData("../out_inv_A.txt", inv_A);

    vector<double> test(SIZE * SIZE);
    Mult(A, inv_A, test, SIZE);
//    WriteData("testA.txt", test);

    vector<double> answer(SIZE);
    Mult_2(inv_A, r, answer, SIZE);
    WriteData("answer_by_gauss.txt", answer);

    // end of first method

    vector<double> inv_L_D(SIZE * SIZE);
    vector<double> inv_L_D__U(SIZE * SIZE);
    vector<double> inv_L_D__f(SIZE);
    Inverse(L_D, inv_L_D, SIZE);
    Mult(inv_L_D, U, inv_L_D__U, SIZE);
    Mult_2(inv_L_D, r, inv_L_D__f, SIZE);

    vector<double> answer_by_iter(SIZE);
    answer_by_iter = Solve_3(A, inv_L_D__U, inv_L_D__f, x_0, r, number_of_iterations, SIZE);
    WriteData("answer_by_iter.txt", answer_by_iter);

//    std::cout << first_norm(A, SIZE) * first_norm(inv_A, SIZE) << std::endl;

//    double max = 0.;
//    double min = 0.;
//    compute_max_and_min_eig(A, SIZE, 20, max, min);
//
//    std::cout << max << " " << min << std::endl;

    vector<double> res_gauss(SIZE);
    vector<double> res_iter(SIZE);
    Mult_2(A, answer, res_gauss, SIZE);
    Mult_2(A, answer_by_iter, res_iter, SIZE);

    for (int i = 0; i < r.size(); i++) {
        res_gauss[i] -= r[i];
        res_iter[i] -= r[i];
    }

    WriteData("res_gauss.txt", res_gauss);
    WriteData("res_iter.txt", res_iter);

    return 0;
}