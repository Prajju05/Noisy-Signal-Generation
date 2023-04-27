#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include "gnuplot-iostream.h"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

const int Fs = 1000; // Sampling frequency
const int N = 1000; // Signal length
const double f1 = 10; // Frequency of sinusoid 1
const double f2 = 20; // Frequency of sinusoid 2
const double f3 = 30; // Frequency of sinusoid 3

int main()
{
  std::vector<double> t(N); // Time vector
std::vector<double> x(N); // Signal vector

for (int n = 0; n < N; n++) {
    t[n] = n / (double)Fs;
    x[n] = sin(2 * M_PI * f1 * t[n]) + 0.5 * sin(2 * M_PI * f2 * t[n]) + 0.3 * sin(2 * M_PI * f3 * t[n]);
}
std::default_random_engine generator;
std::normal_distribution<double> distribution(0, sqrt(pow(10, -20/10.0) / 2)); // -20 dBm RMS noise
std::vector<double> noise(N);

for (int n = 0; n < N; n++) {
    noise[n] = distribution(generator);
}
std::vector<double> y(N);

for (int n = 0; n < N; n++) {
    y[n] = x[n] + noise[n];
}
Eigen::MatrixXd X(N, 1);
Eigen::MatrixXd Y(N, 1);

for (int n = 0; n < N; n++) {
    X(n, 0) = x[n];
    Y(n, 0) = y[n];
}

Eigen::JacobiSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
Eigen::MatrixXd U = svd.matrixU();
Eigen::MatrixXd V = svd.matrixV();
Eigen::VectorXd S = svd.singularValues();
double sigma = S(0) * pow(10, -3); // Choose threshold as 3 orders of magnitude below largest singular value
Eigen::MatrixXd D(N, N);
D.setZero();

for (int i = 0; i < N; i++) {
    if (S(i) > sigma) {
        D(i, i) = 1 / S(i);
    }
}

Eigen::MatrixXd X_tilde = V * D * U.transpose() * Y;
std::vector<double> x_tilde(N);

for (int n = 0; n < N; n++) {
    x_tilde[n] = X_tilde(n, 0);
}
// Plot the signals
Gnuplot gp;
gp << "set multiplot layout 2,2\n";
gp << "set xlabel 'Time (s)'\n";
gp << "set ylabel '

}