    #include <iostream>
#include <vector>
#include <cmath> // for fabs
#include <iomanip> // for std::setprecision
using namespace std;

const double EPSILON = 1e-9; // Threshold to check for near-zero pivots

// Function to perform Gauss-Jordan Elimination
void gaussJordan(vector<vector<double>> A, vector<double>& B) {
    int n = A.size();

    // Augmented matrix
    for (int i = 0; i < n; i++) {
        A[i].push_back(B[i]);
    }

    // Perform the elimination process
    for (int i = 0; i < n; i++) {
        // Partial Pivoting
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(A[i], A[maxRow]);

        // Check for nearly singular matrix
        if (fabs(A[i][i]) < EPSILON) {
            cerr << "Matrix is singular or nearly singular!" << endl;
            return;
        }

        // Normalize the pivot row
        double pivot = A[i][i];
        for (int j = 0; j <= n; j++) {
            A[i][j] /= pivot;
        }

        // Eliminate all other rows
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = 0; j <= n; j++) {
                    A[k][j] -= factor * A[i][j];
                }
            }
        }
    }
    // for (auto l : A) {for(auto k : l) cout << k << " "; cout << endl;}
    // Extract the solution
    for (int i = 0; i < n; i++) {
        B[i] = A[i][n];
    }
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

    gaussJordan(A, B);

    cout << "Solution: " << endl;
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << setprecision(6) << B[i] << endl;
    }

    return 0;
}
