#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

const double EPSILON = 1e-9;
// by chatgpt

// Function to perform Gaussian Elimination
vector<double> gaussianElimination(vector<vector<double>> A, vector<double>& B) {
    int n = A.size();

    // Augmented matrix
    for (int i = 0; i < n; i++) {
        A[i].push_back(B[i]);
    }

    // Forward Elimination
    for (int i = 0; i < n; i++) {
        // Partial Pivoting
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(A[i], A[maxRow]);

        // Check for singular matrix
        if (fabs(A[i][i]) < EPSILON) {
            cerr << "Matrix is singular or nearly singular!" << endl;
            return {};
        }

        // Eliminate entries below the pivot
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j <= n; j++) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }

    // for (auto l : A) {for(auto k : l) cout << k << " "; cout << endl;}

    // Back Substitution
    vector<double> solution(n);
    for (int i = n - 1; i >= 0; i--) {
        solution[i] = A[i][n] / A[i][i];
        for (int j = i - 1; j >= 0; j--) {
            A[j][n] -= A[j][i] * solution[i];
        }
    }

    return solution;
}

int main() {
    int n;
    cout << "Enter the number of variables: ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> B(n);

    cout << "Enter the coefficients of the matrix A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    cout << "Enter the constants of the vector B:" << endl;
    for (int i = 0; i < n; i++) {
        cin >> B[i];
    }

    vector<double> solution = gaussianElimination(A, B);

    if (!solution.empty()) {
        cout << "Solution: " << endl;
        for (int i = 0; i < n; i++) {
            cout << "x" << i + 1 << " = " << solution[i] << endl;
        }
    }

    return 0;
}
