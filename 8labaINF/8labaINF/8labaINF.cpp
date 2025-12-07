#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;

typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;
const double EPSILON = 1e-9;
const int MAX_ITER = 100;
const int N = 3;

void printVector(const Vector& V) {
    for (int i = 0; i < N; i++) {
        cout << "x" << i + 1 << " = " << setw(10) << setprecision(6) << fixed << V[i] << endl;
    }
    cout << endl;
}

void printMatrix(const Matrix& A) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << setw(10) << setprecision(6) << fixed << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void printTVerh(const string& method, double epsilon) {
    cout << "\nМетод " << method << " Eps " << setprecision(5) << epsilon << " Лимит " << MAX_ITER << endl;
    cout << setw(3) << "i"
        << setw(13) << "x1"
        << setw(13) << "x2"
        << setw(13) << "x3" << endl;
}

void printT(int i, const Vector& X) {
    cout << setw(3) << i;
    cout << setprecision(6) << fixed << setw(13) << X[0];
    cout << setprecision(6) << fixed << setw(13) << X[1];
    cout << setprecision(6) << fixed << setw(13) << X[2] << endl;
}

double findDelta(const Vector& X_new, const Vector& X_old) {
    double max_delta = 0.0;
    for (int i = 0; i < N; ++i) {
        max_delta = max(max_delta, abs(X_new[i] - X_old[i]));
    }
    return max_delta;
}

void readMatrix(Matrix& A, Vector& B) {
    cout << "\nВвод Системы 3x3" << endl;

    cout << "Матрица A:" << endl;
    for (int i = 0; i < N; i++) {
        cout << "Строка " << i + 1 << ": ";
        for (int j = 0; j < N; j++) {
            cin >> A[i][j];
        }
    }

    cout << "Вектор B:" << endl;
    for (int i = 0; i < N; i++) {
        cout << "B[" << i + 1 << "]: ";
        cin >> B[i];
    }
}

Matrix invMatrix(Matrix A) {
    int n = A.size();
    Matrix I(n, Vector(n, 0.0));
    for (int i = 0; i < n; i++) I[i][i] = 1.0;

    for (int i = 0; i < n; i++) {
        double p = A[i][i];
        int pRow = i;
        if (abs(p) < EPSILON) {
            for (int k = i + 1; k < n; k++) {
                if (abs(A[k][i]) > abs(p)) {
                    p = A[k][i];
                    pRow = k;
                }
            }
            if (abs(p) < EPSILON) {
                cerr << "Ошибка: Решения нет." << endl;
                return Matrix();
            }
            swap(A[i], A[pRow]);
            swap(I[i], I[pRow]);
            p = A[i][i];
        }

        for (int j = 0; j < n; j++) {
            A[i][j] /= p;
            I[i][j] /= p;
        }

        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = 0; j < n; j++) {
                    A[k][j] -= factor * A[i][j];
                    I[k][j] -= factor * I[i][j];
                }
            }
        }
    }
    return I;
}

Vector umnMatrix(const Matrix& A, const Vector& B) {
    int n = A.size();
    Vector Result(n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Result[i] += A[i][j] * B[j];
        }
    }
    return Result;
}

void resInv(const Matrix& A, const Vector& B) {
    cout << "\nРешение Матричным Методом" << endl;
    Matrix A_inv = invMatrix(A);
    if (A_inv.empty()) return;

    cout << "Обратная Матрица A-1:" << endl;
    printMatrix(A_inv);

    Vector X = umnMatrix(A_inv, B);
    cout << "Вектор решений X:" << endl;
    printVector(X);
}

void Gaus(Matrix A, Vector B) {
    cout << "\nРешение Методом Гаусса" << endl;

    for (int k = 0; k < N; k++) {
        int pRow = k;
        for (int i = k + 1; i < N; i++) {
            if (abs(A[i][k]) > abs(A[pRow][k])) pRow = i;
        }
        swap(A[k], A[pRow]);
        swap(B[k], B[pRow]);

        if (abs(A[k][k]) < EPSILON) {
            cerr << "Ошибка: Решения нет." << endl;
            return;
        }

        for (int i = k + 1; i < N; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < N; j++) {
                A[i][j] -= factor * A[k][j];
            }
            B[i] -= factor * B[k];
        }
    }

    Vector X(N);
    for (int i = N - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < N; j++) {
            sum += A[i][j] * X[j];
        }
        X[i] = (B[i] - sum) / A[i][i];
    }

    cout << "Вектор решений X:" << endl;
    printVector(X);
}

void resIter(const Matrix& A, const Vector& B, double epsilon, bool isSeidel) {
    string method = isSeidel ? "ЗЕЙДЕЛЯ" : "ЯКОБИ";
    printTVerh(method, epsilon);

    Vector X_old(N, 0.0);
    Vector X_new = X_old;
    double max_delta = 0.0;
    int k = 0;

    printT(k, X_old);

    for (k = 1; k <= MAX_ITER; ++k) {

        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j) {
                if (i != j) {
                    if (isSeidel && j < i) {
                        sum += A[i][j] * X_new[j];
                    }
                    else {
                        sum += A[i][j] * X_old[j];
                    }
                }
            }

            X_new[i] = (B[i] - sum) / A[i][i];
        }

        max_delta = findDelta(X_new, X_old);
        printT(k, X_new);

        if (max_delta < epsilon) {
            cout << "Решение найдено на шаге k = " << k << endl;
            return;
        }

        X_old = X_new;
    }

    cout << "Достигнут лимит " << MAX_ITER << " итераций. Метод расходится." << endl;
}

void printIter(const Matrix& A, const Vector& B, bool isSeidel) {
    cout << "\nEps = 0.01" << endl;
    resIter(A, B, 0.01, isSeidel);

    cout << "\nEps = 0.00001" << endl;
    resIter(A, B, 0.00001, isSeidel);
}


int main() {
    setlocale(LC_ALL, "Russian");

    Matrix A(N, Vector(N));
    Vector B(N);

    readMatrix(A, B);

    int choice;
    do {
        cout << "\nТекущая система A:" << endl;
        printMatrix(A);
        cout << "Вектор B:";
        for (double val : B) cout << " " << val;
        cout << "\n\nВыбор метода:" << endl;

        cout << "1. Матричный метод (A-1 * B)" << endl;
        cout << "2. Метод Гаусса" << endl;
        cout << "3. Метод Якоби" << endl;
        cout << "4. Метод Зейделя" << endl;
        cout << "5. Ввести другую систему" << endl;
        cout << "0. Выход" << endl;
        cin >> choice;

        if (cin.fail()) {
            cin.clear(); 
            cin.ignore(10000, '\n');
            choice = -1; 
        }

        switch (choice) {
        case 1:
            resInv(A, B);
            break;
        case 2:
            Gaus(A, B);
            break;
        case 3:
            printIter(A, B, false);
            break;
        case 4:
            printIter(A, B, true);
            break;
        case 5:
            readMatrix(A, B);
            break;
        case 0:
            cout << "Программа завершена." << endl;
            break;
        default:
            cout << "Неверный выбор." << endl;
        }
    } while (choice != 0);

    return 0;
}