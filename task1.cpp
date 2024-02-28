#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char* argv[])
{
    for (int n = 10; n < 100000001; n *= 10) {
        double h = 1.0 / n;
        double* ans = new double[n + 1];
        double* alpha = new double[n];
        double* beta = new double[n];
        double* f = new double[n];
        // -u'' = sin(x)
        // u(0) = 0
        // u(1) = sin(1) ~ 0.84171
        // u(x) = sin(x) - exact solution
        ans[0] = 0;       
        ans[n] = sin(1); 
        f[0] = 0;
        f[n] = sin(1);
        // A_i = C_i = 1, B_i = -2, for all i
        alpha[0] = -1/2; // alpha_0 = -C1/B1
        beta[0] = 0;     // beta_0  =  F1/B1
        // alpha_i+1 = -C_i / (A_i * alpha_i + B_i)
        // beta_i+1 = (F_i - A_i * beta_i) / (A_i * alpha_i + B_i)
        int A = -1;
        int C = -1;
        int B = 2;
        for (int i = 1; i < n; i++) {
            f[i] = sin(i * h) * h * h;
            alpha[i] = -C / (A * alpha[i - 1] + B);
            beta[i] = (f[i - 1] - A * beta[i - 1]) / (A * alpha[i - 1] + B);
        }
        f[n - 1] += sin(1);
        // ans_n-1 = (F_n-1 - A_n-1 * beta_n-1) / (B_n-1 + A_n-1 * alpha_n-1)
        ans[n - 1] = (f[n - 1] - A * beta[n - 1]) / (B + A * alpha[n - 1]);
        for (int i = n - 2; i > 0; i--) {
            ans[i] = alpha[i + 1] * ans[i + 1] + beta[i + 1];
        }
        double C_error = 0;
        double L2_error = 0;
        double exact = 0;
        for (int i = 1; i < n; i++) {
            exact = sin(i * h);
            if (C_error < abs(ans[i] - exact)) {
                C_error = abs(ans[i] - exact);
            }
            L2_error += (ans[i] - exact) * (ans[i] - exact);
        }
        cout << n << endl;
        cout << "L2_error = " << sqrt(L2_error * h) << endl;
        cout << "C_error = " << C_error << endl;
        cout << endl;
        delete[] ans;
        delete[] alpha;
        delete[] beta;
        delete[] f;
    }
        return 0;
}
