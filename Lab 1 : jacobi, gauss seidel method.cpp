#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;

vector<double> coeff[500];
double C[500];
const double eps = 1e-6;

int N = 3;
bool input = false; // make it true for taking custom input for the variable size

double getvalue(int eqnNo, vector<double>& var) { // equation no, current variables value
    // x = (k - (by + cz)) / N
    double ret = C[eqnNo];
    for (int i = 0; i < N; i++) {
        if(i != eqnNo)
            ret -= (coeff[eqnNo][i] * var[i]);
    }
    return ret / coeff[eqnNo][eqnNo];
}

void jacobi() {
    vector<double> prev(N, 0.0), cur(N);

    // 1st iteration
    for (int i = 0; i < N; i++) {
        cur[i] = getvalue(i, prev);
        // cout << cur[i] << " \n"[i == N - 1];
    }
    int iterations = 10;
    vector<double> d(N, 0.0); // like dx, dy, dz
    double eavg;
    while (--iterations) {
        eavg = 0.0;
        swap(prev, cur);
        for (int i = 0; i < N; i++) {
            cur[i] = getvalue(i, prev);
            d[i] = fabs(cur[i] - prev[i]);
            // cout << cur[i] << " \n"[i == N - 1];
            // printf("%lf ", d[i]);
            // if(i == N -1) puts("");
            eavg += d[i];
        }
        eavg = eavg / N;
        printf("%lf\n", eavg);
    }
}

void gauss_seidel() {
    vector<double> cur(N, 0), prev(N, 0);
    // 1st iteration
    for (int i = 0; i < N; i++) {
        cur[i] = getvalue(i, cur);
        // cout << cur[i] << " \n"[i == N - 1];
    }

    int iterations = 10;

    double eavg = 0;
    vector<double> d(N);
    while (--iterations) {
        eavg = 0;
        prev = cur;
        for (int i = 0; i < N; i++) {
            cur[i] = getvalue(i, cur);
            d[i] = fabs(cur[i] - prev[i]);
            eavg += d[i];
            // cout << cur[i] << " \n"[i == N - 1];
            // printf("%lf ", d[i]);
            if(i == N -1) puts("");
        }
        eavg /= N;
        printf("%lf ", eavg);
    }
}

int main() {
    if (input) {
        cout << "Number of variables: ";
        cin >> N;
    }

    // cout << "Next N lines give input in the form\n";
    // cout << "N + 1 inputs for the next N lines\n";
    // cout << "a1 a2 ... aN = c\n";
    for (int i = 0; i < N; i++) {
        double x;
        for (int j = 0; j < N; j++) {
            cin >> x;
            coeff[i].push_back(x);
        }
        cin >> C[i];
    }

    jacobi();
    // gauss_seidel();

    return 0;
}
