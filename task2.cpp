#include "inmost.h"
#include <stdio.h>


using namespace INMOST;
using namespace std;


int main(int argc, char *argv[])
{
    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector sol;
    int N = 100;

    A.SetInterval(0, (N - 1) * (N - 1));
    b.SetInterval(0, (N - 1) * (N - 1));
    sol.SetInterval(0, (N - 1) * (N - 1));

    Solver S(Solver::INNER_ILU2);
    S.SetParameter("absolute_tolerance", "1e-10");
    S.SetParameter("relative_tolerance", "1e-6");

    for(int i = 0; i < (N - 1) * (N - 1); i++) {
        A[i][i] = 4;
        b[i] = 50.0 * sin(5.0 * (double) (i % (N - 1) + 1) / N) * sin(5.0 * (double) (i / (N - 1) + 1) / N) / N / N;           

        if (i % (N - 1) + 1 < (N - 1)) {
            A[i][i + 1] = -1;
        } else {
            b[i] += sin(5.0) * sin(5.0 * (double) (i / (N - 1) + 1) / N);
        }

        if (i % (N - 1) > 0) {
            A[i][i - 1] = -1;
        }

        if (i + N - 1 < (N - 1) * (N - 1)) {
            A[i][i + N - 1] = -1;
        }
        else {
            b[i] += sin(5.0) * sin(5.0 * (double) (i % (N - 1) + 1) / N);
        }

        if (i - N > 0) {
            A[i][i - N + 1] = -1;
        }
    }

    S.SetMatrix(A);

    bool solved = S.Solve(b, sol);

    std::cout << "N = " << N << std::endl;
    double norm = 0;
    for (int i = 0; i < (N - 1) * (N - 1); i++) {
        double temp = abs(sin(5.0 * (double) (i % (N - 1) + 1) / N) * sin(5.0 * (double) (i / (N - 1) + 1) / N) - sol[i]);
        if (temp > norm) {
            norm = temp;
        }
    }
    cout << "C-norm = " << norm << endl;

    norm = 0;
    for (int i = 0; i < (N - 1) * (N - 1); i++) {
        norm += pow((sin(5.0 * (double)(i % (N - 1) + 1) / N) * sin(5.0 * (double)(i / (N - 1) + 1) / N) - sol[i]), 2);
    }
    norm = sqrt(norm) / (N - 1) / (N - 1);
    cout << "L2-norm = " << norm << endl;
    cout << "num.iters: " << S.Iterations() << endl;
    cout << "prec.time: " << S.PreconditionerTime() << endl;
    cout << "iter.time: " << S.IterationsTime() << endl;
    if(!solved){
        cout << "Linear solver failure!" << endl;
        cout << "Reason: " << S.ReturnReason() << endl;
    }
	return 0;
}