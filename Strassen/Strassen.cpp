#include <iostream>
#include <complex>

const int N = 4;

// Оголошення структури для комплексного числа
struct Complex {
    double real;
    double imag;
};

// Функція додавання комплексних чисел
Complex add(Complex a, Complex b) {
    Complex result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

// Функція віднімання комплексних чисел
Complex sub(Complex a, Complex b) {
    Complex result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}

// Функція множення комплексних чисел
Complex multiply(Complex a, Complex b) {
    Complex result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}


// Класичний алгоритм множення матриць
void СlassicMultiply(Complex A[][N], Complex B[][N], Complex C[][N]) {
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            C[i][j].real = 0;
            C[i][j].imag = 0;
            for (int t = 0; t < 2; t++) {
                C[i][j] = add(C[i][j], multiply(A[i][t], B[t][j]));
            }
        }
    }
}

// Додавання двох матриць
void AddMatrix(int n, Complex X[][N], Complex Y[][N], Complex Z[][N]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Z[i][j] = add(X[i][j], Y[i][j]);
        }
    }
}

// Віднімання двох матриць
void SubMatrix(int n, Complex X[][N], Complex Y[][N], Complex Z[][N]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Z[i][j] = sub(X[i][j], Y[i][j]);
        }
    }
}


// Власне, алгоритм Штрассена
void Strassen(int n, Complex A[][N], Complex B[][N], Complex C[][N]) {
    Complex A11[N][N], A12[N][N], A21[N][N], A22[N][N];
    Complex B11[N][N], B12[N][N], B21[N][N], B22[N][N];
    Complex C11[N][N], C12[N][N], C21[N][N], C22[N][N];
    Complex M1[N][N], M2[N][N], M3[N][N], M4[N][N], M5[N][N], M6[N][N], M7[N][N];
    Complex AA[N][N], BB[N][N];

    if (n == 2) {
        СlassicMultiply(A, B, C);
    }
    else {
        //ініціалізація допоміжних матриць
        for (int i = 0; i < n / 2; i++) {
            for (int j = 0; j < n / 2; j++) {
                A11[i][j] = A[i][j];
                A12[i][j] = A[i][j + n / 2];
                A21[i][j] = A[i + n / 2][j];
                A22[i][j] = A[i + n / 2][j + n / 2];

                B11[i][j] = B[i][j];
                B12[i][j] = B[i][j + n / 2];
                B21[i][j] = B[i + n / 2][j];
                B22[i][j] = B[i + n / 2][j + n / 2];
            }
        }

        //обчислення M1 = (A11 + A22) × (B11 + B22)
        AddMatrix(n / 2, A11, A22, AA);
        AddMatrix(n / 2, B11, B22, BB);
        Strassen(n / 2, AA, BB, M1);

        //обчислення M2 = (A21 + A22) × B11
        AddMatrix(n / 2, A21, A22, AA);
        Strassen(n / 2, AA, B11, M2);

        //обчислення M3 = A11 × (B12 - B22)
        SubMatrix(n / 2, B12, B22, BB);
        Strassen(n / 2, A11, BB, M3);

        //обчислення M4 = A22 × (B21 - B11)
        SubMatrix(n / 2, B21, B11, BB);
        Strassen(n / 2, A22, BB, M4);

        //обчислення M5 = (A11 + A12) × B22
        AddMatrix(n / 2, A11, A12, AA);
        Strassen(n / 2, AA, B22, M5);

        //обчислення M6 = (A21 - A11) × (B11 + B12)
        SubMatrix(n / 2, A21, A11, AA);
        AddMatrix(n / 2, B11, B12, BB);
        Strassen(n / 2, AA, BB, M6);

        //обчислення M7 = (A12 - A22) × (B21 + B22)
        SubMatrix(n / 2, A12, A22, AA);
        AddMatrix(n / 2, B21, B22, BB);
        Strassen(n / 2, AA, BB, M7);

        //обчислення C11 = M1 + M4 - M5 + M7
        AddMatrix(n / 2, M1, M4, AA);
        SubMatrix(n / 2, M7, M5, BB);
        AddMatrix(n / 2, AA, BB, C11);

        //обчислення C12 = M3 + M5
        AddMatrix(n / 2, M3, M5, C12);

        //обчислення C21 = M2 + M4
        AddMatrix(n / 2, M2, M4, C21);

        //обчислення C22 = M1 - M2 + M3 + M6
        SubMatrix(n / 2, M1, M2, AA);
        AddMatrix(n / 2, M3, M6, BB);
        AddMatrix(n / 2, AA, BB, C22);

        //злиття в C[][N]
        for (int i = 0; i < n / 2; i++) {
            for (int j = 0; j < n / 2; j++) {
                C[i][j] = C11[i][j];
                C[i][j + n / 2] = C12[i][j];
                C[i + n / 2][j] = C21[i][j];
                C[i + n / 2][j + n / 2] = C22[i][j];
            }
        }
    }
}


// Функція виведення матриці комплексних чисел
void СoutMatrix(int n, Complex C[][N]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << C[i][j].real << " + " << C[i][j].imag << "i" << "\t";
        }
        std::cout << std::endl;
    }

    std::cout << "--------------------------------------------------------" << std::endl;
}

int main() {
    // Цілі додатні числа
    Complex A1[N][N] = {
        { {7, 0}, {3, 0}, {5, 0}, {8, 0} },
        { {1, 0}, {3, 0}, {5, 0}, {6, 0} },
        { {9, 0}, {9, 0}, {0, 0}, {2, 0} },
        { {4, 0}, {6, 0}, {8, 0}, {3, 0} }
    };

    Complex B1[N][N] = {
        { {1, 0}, {7, 0}, {5, 0}, {3, 0} },
        { {9, 0}, {5, 0}, {4, 0}, {7, 0} },
        { {8, 0}, {9, 0}, {1, 0}, {3, 0} },
        { {5, 0}, {2, 0}, {0, 0}, {1, 0} }
    };

    Complex C1[N][N];
    Strassen(N, A1, B1, C1);
    СoutMatrix(N, C1);

    // Уявні числа
    Complex A2[N][N] = {
        { {0, 7}, {0, 3}, {0, 5}, {0, 8} },
        { {0, 1}, {0, 3}, {0, 5}, {0, 6} },
        { {0, 9}, {0, 9}, {0, 0}, {0, 2} },
        { {0, 4}, {0, 6}, {0, 8}, {0, 3} }
    };

    Complex B2[N][N] = {
        { {0, 1}, {0, 7}, {0, 5}, {0, 3} },
        { {0, 9}, {0, 5}, {0, 4}, {0, 7} },
        { {0, 8}, {0, 9}, {0, 1}, {0, 3} },
        { {0, 5}, {0, 2}, {0, 0}, {0, 1} }
    };

    Complex C2[N][N];
    Strassen(N, A2, B2, C2);
    СoutMatrix(N, C2);

    // Комплексні числа
    Complex A3[N][N] = {
        { {7, 7}, {3, 3}, {5, 5}, {8, 8} },
        { {1, 1}, {3, 3}, {5, 5}, {6, 6} },
        { {9, 9}, {9, 9}, {0, 0}, {2, 2} },
        { {4, 4}, {6, 6}, {8, 8}, {3, 3} }
    };

    Complex B3[N][N] = {
        { {1, 1}, {7, 7}, {5, 5}, {3, 3} },
        { {9, 9}, {5, 5}, {4, 4}, {7, 7} },
        { {8, 8}, {9, 9}, {1, 1}, {3, 3} },
        { {5, 5}, {2, 2}, {0, 0}, {1, 1} }
    };

    Complex C3[N][N];
    Strassen(N, A3, B3, C3);
    СoutMatrix(N, C3);

    // Комплексні числа з різними знаками + неквадратна матриця
    Complex A4[N][N] = {
        { {-7.3, 0}, {3, 0}},
        { {0, -3.4}, {3, 0}},
        { {2, 0}, {3, 0}},
        { {0, 0}, {0, -1}},

    };

    Complex B4[N][N] = {
        { {1, 0}, {7, 0}, {5, 0}, {3, 0} },
        { {9, 0}, {5, 0}, {4, 0}, {7, 0} },
        { {8, 0}, {9, 0}, {1, 0}, {3, 0} },
        { {5, 0}, {2, 0}, {0, 0}, {1, 0} }
    };

    Complex C4[N][N];
    Strassen(N, A4, B4, C4);
    СoutMatrix(N, C4);

    return 0;
}
