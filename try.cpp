#include <bits/stdc++.h>

using namespace std;

typedef complex<double> cd;

vector<cd> fft(vector<cd> a, bool invert = false) {
    int n = a.size();

    if (n == 1) {
        return a;
    }

    vector<cd> w(n);
    for (int i = 0; i < n; ++i) {
        // Sign is negative for forward FFT, positive for inverse IFFT
        double alpha = 2 * M_PI * i / n * (invert ? 1 : -1);
        w[i] = cd(cos(alpha), sin(alpha));
	}

    vector<cd> ae(n / 2);
    vector<cd> ao(n / 2);

    for (int i = 0; i < n / 2; i++) {
        ae[i] = a[2 * i];
        ao[i] = a[2 * i + 1];
    }

    vector<cd> ye = fft(ae, invert);
    vector<cd> yo = fft(ao, invert);

    vector<cd> y(n);
    for (int k = 0; k < n / 2; ++k) {
        y[k] = ye[k] + w[k] * yo[k];
        y[k + n / 2] = ye[k] - w[k] * yo[k];

	}

    return y;
}

vector<double> multiply_polynomials(vector<double>& poly1, vector<double>& poly2) {
    int n = 1;
    while (n < poly1.size() + poly2.size()) {
        n <<= 1; 
    }

    vector<cd> a(n, 0), b(n, 0);
    for (size_t i = 0; i < poly1.size(); i++) a[i] = poly1[i];
    for (size_t i = 0; i < poly2.size(); i++) b[i] = poly2[i];

    vector<cd> fa = fft(a, false);
    vector<cd> fb = fft(b, false);

    vector<cd> fc(n);
    for (int i = 0; i < n; i++) {
        fc[i] = fa[i] * fb[i];
    }

    vector<cd> c = fft(fc, true);

    vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = round(c[i].real() / n);
    }

    while (result.size() > 1 && result.back() == 0) {
        result.pop_back();
    }

    return result;
}

int main() {
    vector<double> poly1 = {1.0, 2.0, 3.0};
    vector<double> poly2 = {4.0, 5.0};

    vector<double> product = multiply_polynomials(poly1, poly2);

    cout << "Resulting Polynomial Coefficients: ";
    for (double coef : product) {
        cout << coef << " ";
    }
    cout << endl;

    return 0;
}
