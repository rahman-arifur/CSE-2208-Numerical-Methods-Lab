#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

// Md. Arifur Rahman 
// 2107062
const double tolerance = 1e-5;

template <typename T>
T power(T b, int p) {
    T res = 1;
    for(;p > 0; p >>= 1, b *= b)
        if(p & 1) res *= b;
    return res;
}

double getval(double x, vector<double>& a) {
    double res = 0;
    for (int i = 0; i < size(a); ++i) {
        res += power(x, i) * a[i];
    }
    return res;
}

double fprimex(double x, vector<double>& a) {
    double res = 0;
    for (int i = 1; i < size(a); ++i) {
        res += power(x, i - 1) * a[i] * i;
    }
    return res;
}

void bisection(vector<double>& a) {
    int n = size(a) - 1;
    double val = power(a[n-1] / a[n], 2) - 2 * a[n - 2] / a[n]; 
    val = fabs(val);
    double left = sqrt(val), right = -left, c;
    
    int it = 0;
    while (fabs(left - right) > tolerance and it < 100) {
        c = (left + right) / 2;
        ++it;
        double fc = getval(c, a), fa = getval(left, a);
        if (fc == 0) break;
        else if (fa * fc < 0) right = c;
        else left = c;
    }
    printf("%.6lf\n", c);
}

double falsePosition(double x1, double x2, vector<double>& a) {
    return (x1 * getval(x2, a) - x2 * getval(x1, a)) / (getval(x2, a) - getval(x1, a));
}
void falsePosition(vector<double>& a) {
    int n = size(a) - 1;
    double left, right, c;
    cout << "Beginning a, b:\n";
    cin >> left >> right;
    int it = 0;
    if (getval(left, a) * getval(right, a) >= 0) {
        cout << "Enter 2 valid values such that f(a) * f(b) < 0\n";
        return;
    }
    while(it < 100 and getval(c, a) <= tolerance) {
        c = falsePosition(left, right, a);
        if (getval(c, a) == 0) break;
        if (getval(left, a) * getval(c, a) < 0) right = c;
        else left = c;
    }
    printf("%.6lf\n", c);
}

void newtonRaphson(vector<double>& a) {
    double x0, x1, ini;
    cout << "Beginning value: ";
    cin >> x0;
    ini = x0;
    int it = 0;
    while (it < 100) {
        double k = fprimex(x0, a);
        if (k == 0) {
            cout << "Can't solve with initial value = " << ini << "\n";
            return;
        }
        x1 = x0 - getval(x0, a) / k;
        if (abs(x1 - x0) < tolerance)
            break;
        ++it;
        x0 = x1;
    }
    printf("%.6lf\n", x1);
}

void secant(vector<double>& a) {
    double x1, x2, x3;
    cout << "Two initial guess range: ";
    cin >> x1 >> x2;

    if (getval(x1, a) * getval(x2, a) >= 0) {
        cout << "Secant method cannot be applied for " << x1 << " and " << x2 << '\n';
        return;
    }

    int it = 0;
    while (it < 100) {
        double fx1 = getval(x1, a), fx2 = getval(x2, a);
        x3 = fx2 * (x2 - x1) / (fx2 - fx1);
        x3 = x2 - x3;
        if (fabs(getval(x3, a)) <= tolerance) break;
        x1 = x2;
        x2 = x3;
        ++it;
    }
    printf("%.6lf\n", x3);
}

int main() 
{
    int n = 2; // maximum power of x
    cout << "Max power of x: ";
    cin >> n;
    vector<double> a(n + 1);
    cout << "input in the form anx^n + a(n-1)x^(n-1) + .. + a1x + a0 = 0\n";
    for (int i = n; i >= 0; --i) {
        cin >> a[i];
    }
    
    // bisection(a);
    // falsePosition(a);
    // newtonRaphson(a);
    // secant(a);
return 0;
}