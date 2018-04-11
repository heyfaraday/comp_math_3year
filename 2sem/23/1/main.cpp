#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

double f(double x_n, double z_n, double epsilon)
{
    double mean = -x_n + epsilon * z_n * (1. - x_n * x_n);
    return mean;
}

double beta(double p)
{
    return pow((0.25 / 0.38), 1. / (p + 1.));
}

void write_result(const char *name_x, const char *name_z, const std::vector<std::vector<double>> &result, const std::vector<double> &grid)
{

    FILE *file_x = fopen(name_x, "w");
    FILE *file_z = fopen(name_z, "w");

    if (file_x == 0)
    {
        printf("Cannot create file: %s\n", name_x);
        return;
    }

    if (file_z == 0)
    {
        printf("Cannot create file: %s\n", name_z);
        return;
    }

    for (int i = 0; i < result.size(); ++i)
    {
        fprintf(file_x, "%.6e  %.6e\n", grid.at(i), result.at(i).at(0));
        fprintf(file_z, "%.6e  %.6e\n", grid.at(i), result.at(i).at(1));
    }

    fclose(file_x);
    fclose(file_z);
}

std::vector<double> make_runge_kutta_4_step(const std::vector<double> &x_z_n, double h, double epsilon)
{

    std::vector<double> x_z_n_1 = {0., 0.};

    double k_0 = x_z_n.at(1);
    double l_0 = f(x_z_n.at(0), x_z_n.at(1), epsilon);

    double k_1 = x_z_n.at(1) + 0.5 * l_0 * h;
    double l_1 = f(x_z_n.at(0) + 0.5 * h * k_0, x_z_n.at(1) + 0.5 * h * l_0, epsilon);

    double k_2 = x_z_n.at(1) + 0.5 * l_1 * h;
    double l_2 = f(x_z_n.at(0) + 0.5 * h * k_1, x_z_n.at(1) + 0.5 * h * l_1, epsilon);

    double k_3 = x_z_n.at(1) + l_2 * h;
    double l_3 = f(x_z_n.at(0) + h * k_2, x_z_n.at(1) + l_2 * h, epsilon);

    x_z_n_1.at(0) = x_z_n.at(0) + h * (1. / 6.) * (k_0 + 2. * k_1 + 2. * k_2 + k_3);
    x_z_n_1.at(1) = x_z_n.at(1) + h * (1. / 6.) * (l_0 + 2. * l_1 + 2. * l_2 + l_3);

    return x_z_n_1;
}

std::vector<double> make_runge_kutta_1_step(const std::vector<double> &x_z_n, double h, double epsilon)
{

    std::vector<double> x_z_n_1 = {0., 0.};

    double k_0 = x_z_n.at(1);
    double l_0 = f(x_z_n.at(0), x_z_n.at(1), epsilon);

    x_z_n_1.at(0) = x_z_n.at(0) + h * k_0;
    x_z_n_1.at(1) = x_z_n.at(1) + h * l_0;

    return x_z_n_1;
}

void runge_kutta_1_iter_with_precision(std::vector<std::vector<double>> &result,
                                       std::vector<double> &grid,
                                       const std::vector<double> &x_z_0,
                                       double h_0,
                                       double T,
                                       double epsilon,
                                       double t_up)
{

    result.push_back(x_z_0);
    grid.push_back(0.);

    double error = 0.;
    double h = h_0;

    double current_t = 0.;
    while (current_t < t_up)
    {

        std::vector<double> x_z_next = make_runge_kutta_1_step(result.back(), h, epsilon);
        std::vector<double> x_z_next_half = make_runge_kutta_1_step(result.back(), h / 2., epsilon);
        std::vector<double> x_z_next_2n = make_runge_kutta_1_step(x_z_next_half, h / 2., epsilon);

        error = abs((x_z_next.at(0) - x_z_next_2n.at(0)) / (pow(2., 2. - 1.) - 1.));

        if (error < T)
        {
            result.push_back(x_z_next);
            grid.push_back(current_t + h);

            current_t += h;
        }
        else
        {
            h = beta(2.) * h * pow(T / error, 1. / (2. + 1.));
        }
    }
}

void runge_kutta_4_iter_with_precision(std::vector<std::vector<double>> &result,
                                       std::vector<double> &grid,
                                       const std::vector<double> &x_z_0,
                                       double h_0,
                                       double T,
                                       double epsilon,
                                       double t_up)
{

    result.push_back(x_z_0);
    grid.push_back(0.);

    double error = 0.;
    double h = h_0;

    double current_t = 0.;
    while (current_t < t_up)
    {

        std::vector<double> x_z_next = make_runge_kutta_4_step(result.back(), h, epsilon);
        std::vector<double> x_z_next_half = make_runge_kutta_4_step(result.back(), h / 2., epsilon);
        std::vector<double> x_z_next_2n = make_runge_kutta_4_step(x_z_next_half, h / 2., epsilon);

        error = abs((x_z_next.at(0) - x_z_next_2n.at(0)) / (pow(2., 5. - 1.) - 1.));

        if (error < T)
        {
            result.push_back(x_z_next);
            grid.push_back(current_t + h);

            current_t += h;
        }
        else
        {
            h = beta(5.) * h * pow(T / error, 1. / (5. + 1.));
        }
    }
}

int main()
{

    double t_up = 10.;

    double N = 10;
    double h_0 = 1 / N;
    double T = 0.000001;

    double epsilon = 1.;

    std::vector<double> x_z_0 = {1, 1};

    std::vector<std::vector<double>> result_1;
    std::vector<double> grid_1;
    runge_kutta_1_iter_with_precision(result_1, grid_1, x_z_0, h_0, T, epsilon, t_up);
    write_result("../result_1_x", "../result_1_z", result_1, grid_1);

    std::cout << grid_1.size() << std::endl;

    std::vector<std::vector<double>> result_4;
    std::vector<double> grid_4;

    {
        unsigned int start_time = clock();

        runge_kutta_4_iter_with_precision(result_4, grid_4, x_z_0, h_0, T, epsilon, t_up);

        unsigned int end_time = clock();
        unsigned int search_time = end_time - start_time;

        std::cout << "time: " << search_time << std::endl;
    }

    write_result("../result_4_x", "../result_4_z", result_4, grid_4);

    std::cout << grid_4.size() << std::endl;

    return 0;
}