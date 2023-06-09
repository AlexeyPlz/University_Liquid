#include <iostream>
#include <locale.h>
#include <vector>
#include <math.h>
#include <fstream>
#include <chrono>

using namespace std;
using namespace std::chrono;


// Класс Cужение сосуда
class VesselConstriction {
private:
    double l1, r;

public:
    VesselConstriction() {};
    VesselConstriction(double l1_, double r_) { this->l1 = l1_; this->r = r_; };

    double get_l1() { return l1; }
    double get_r()  { return r; }

    ~VesselConstriction() { cout << "Memory class VesselConstriction cleared." << endl; }
};

// Класс Клапан сосуда
class VesselValve {
private:
    double l2, l3, Rk;

public:
    VesselValve() {};
    VesselValve(double l2_, double l3_, double Rk_) { this->l2 = l2_; this->l3 = l3_; this->Rk = Rk_; };

    double get_l2() { return l2; }
    double get_l3() { return l3; }
    double get_Rk() { return Rk; }

    ~VesselValve() { cout << "Memory class VesselValve cleared." << endl; }
};

// Класс Сетка
class Grid {
private:
    double L, R, h1, h2, tau, q, Re, eps;
    int N1, N2;
    vector<double> x, y;
    vector<vector<double>> c1, c2;
    vector<vector<double>> psi_n, psi_np05, psi_np1;
    vector<vector<double>> omega_n, omega_np05, omega_np1;
    VesselConstriction vc;
    VesselValve vv;

public:
    Grid() {};
    Grid(double L_, double R_, int N1_, int N2_, double q_, double Re_, double eps_, VesselConstriction vc_, VesselValve vv_) {
        this->L = L_; this->R = R_;
        this->N1 = N1_; this->N2 = N2_;
        this->h1 = L_ / N1_; this->h2 = R_ / N2_; 
        this->tau = min(h1 * h1, h2 * h2);
        this->q = q_; this->Re = Re_; this->eps = eps_;
        this->vc = vc_; this->vv = vv_;

        vector<double> x(N1_ + 1), y(N2_ + 1);
        vector<vector<double>> c1(N1_ + 1, vector<double>(N2_ + 1)), c2(N1_ + 1, vector<double>(N2_ + 1));
        vector<vector<double>> psi(N1_ + 1, vector<double>(N2_ + 1)), omega(N1_ + 1, vector<double>(N2_ + 1));

        this->x = x; this->y = y;
        this->c1 = c1; this->c2 = c2;
        this->psi_n = psi, this->psi_np05 = psi, this->psi_np1 = psi;
        this->omega_n = omega, this->omega_np05 = omega, this->omega_np1 = omega;

        for (int i = 0; i <= N1; i++) x[i] = i * h1;
        for (int j = 0; j <= N2; j++) y[j] = j * h2;
        this->x = x; this->y = y;

        for (int i = 0; i <= N1; i++) {
            for (int j = 0; j <= N2; j++) {
                if ((x[i] - vc.get_l1()) * (x[i] - vc.get_l1()) + (y[j] - R) * (y[j] - R) <= vc.get_r() * vc.get_r()) this->c1[i][j] = 1 / eps;
                else this->c1[i][j] = 0;
                if (x[i] >= vv.get_l2() && x[i] <= vv.get_l3() && y[j] <= vv.get_Rk()) this->c2[i][j] = 1 / eps;
                else this->c2[i][j] = 0;
            }
        }

        for (int i = 1; i < N1; i++) {
            for (int j = 1; j < N2; j++) {
                this->psi_n[i][j] = 0; this->omega_n[i][j] = 0;
            }
        }

        for (int j = 1; j < N2; j++) {
            double Pj = psi_P(j), Uj = dU(j);
            this->psi_n[0][j] = Pj; this->psi_n[N1][j] = Pj;
            this->psi_np1[0][j] = Pj; this->psi_np1[N1][j] = Pj;
            this->omega_n[0][j] = 2 / h1 * h1 * Pj - Uj; this->omega_n[N1][j] = 2 / h1 * h1 * Pj - Uj;
        }

        this->psi_n[0][0] = 0, this->psi_n[N1][0] = 0;
        this->psi_np05[0][0] = 0, this->psi_np05[N1][0] = 0;
        this->psi_np1[0][0] = 0, this->psi_np1[N1][0] = 0;
        this->psi_n[0][N2] = q / 2; this->psi_n[N1][N2] = q / 2;
        this->psi_np05[0][N2] = q / 2; this->psi_np05[N1][N2] = q / 2;
        this->psi_np1[0][N2] = q / 2; this->psi_np1[N1][N2] = q / 2;
        this->omega_n[0][0] = -dU(0); this->omega_n[0][N2] = -dU(N2); this->omega_n[N1][0] = -dU(0); this->omega_n[N1][N2] = -dU(N2);
        this->omega_np1[0][0] = -dU(0); this->omega_np1[0][N2] = -dU(N2); this->omega_np1[N1][0] = -dU(0); this->omega_np1[N1][N2] = -dU(N2);

        for (int i = 1; i < N1; i++) {
            this->psi_n[i][0] = 0; this->psi_n[i][N2] = q / 2;
            this->psi_np05[i][0] = 0; this->psi_np05[i][N2] = q / 2;
            this->psi_np1[i][0] = 0; this->psi_np1[i][N2] = q / 2;
            this->omega_n[i][0] = 0; this->omega_n[i][N2] = 2 / h2 * h2 * q / 2;
            this->omega_np1[i][0] = 0; this->omega_np1[i][N2] = 2 / h2 * h2 * q / 2;
        }
    }

    double U(int j) { return 3 * q * (R * R - y[j] * y[j]) / (4 * R * R * R); }
    double dU(int j) { return -3 * q * y[j] / (2 * R * R * R); }
    double psi_P(int j) { return 3 * q * y[j] * (R * R - y[j] * y[j] / 3) / (4 * R * R * R); }

    double get_L() { return L; }
    double get_R() { return R; }
    double get_N1() { return N1; }
    double get_N2() { return N2; }
    double get_h1() { return h1; }
    double get_h2() { return h2; }
    double get_q() { return q; }
    double get_Re() { return Re; }
    double get_eps() { return eps; }
    double get_tau() { return tau; }
    double get_x(int i) { return x[i]; }
    double get_y(int j) { return y[j]; }
    double get_c1(int i, int j) { return c1[i][j]; }
    double get_c2(int i, int j) { return c2[i][j]; }
    double get_psi_n(int i, int j) { return psi_n[i][j]; }
    double get_omega_n(int i, int j) { return omega_n[i][j]; }
    double get_psi_np05(int i, int j) { return psi_np05[i][j]; }
    double get_omega_np05(int i, int j) { return omega_np05[i][j]; }
    double get_psi_np1(int i, int j) { return psi_np1[i][j]; }
    double get_omega_np1(int i, int j) { return omega_np1[i][j]; }

    void set_psi_n(int i, int j, double value) { this->psi_n[i][j] = value; }
    void set_psi_np05(int i, int j, double value) { this->psi_np05[i][j] = value; }
    void set_psi_np1(int i, int j, double value) { this->psi_np1[i][j] = value; }
    void set_omega_n(int i, int j, double value) { this->omega_n[i][j] = value; }
    void set_omega_np05(int i, int j, double value) { this->omega_np05[i][j] = value; }
    void set_omega_np1(int i, int j, double value) { this->omega_np1[i][j] = value; }

    void saveInFile(int it, double norma_max, double delta, double time, VesselConstriction vc, VesselValve vv) {
        ofstream file_result("result.txt");
        file_result << "Последовательный алгоритм: " << endl;
        file_result << "1) Число итераций = " << it << endl;
        file_result << "2) eps = " << eps << endl;
        file_result << "3) Норма = " << norma_max << ",   delta = " << delta << endl;
        file_result << "4) Время выполнения программы = " << time << " s" << endl;

        ofstream file_data("data.txt");
        file_data << L << " " << R << endl;
        file_data << N1 + 1 << " " << 2 * N2 + 1 << endl;
        file_data << vc.get_l1() << " " << vc.get_r() << " " << vv.get_l2() << " " << vv.get_l3() << " " << vv.get_Rk() << endl;

        for (int i = 0; i <= N1; i++) file_data << x[i] << " ";
        file_data << endl;
        for (int j = 0; j <= 2 * N2; j++) file_data << -R + j * h2 << " ";
        file_data << endl;

        vector<vector<double>> psi(N1 + 1, vector<double>(2 * N2 + 1));
        vector<vector<double>> omega(N1 + 1, vector<double>(2 * N2 + 1));

        for (int i = 0; i <= N1; i++) {
            for (int j = 0; j < N2; j++) {
                psi[i][j] = -psi_np1[i][N2 - j];
                omega[i][j] = -omega_np1[i][N2 - j];
            }
            for (int j = N2; j <= 2 * N2; j++) {
                psi[i][j] = psi_np1[i][j - N2];
                omega[i][j] = omega_np1[i][j - N2];
            }
        }

        ofstream file_psi("psi.txt");
        for (int j = 0; j <= 2 * N2; j++) {
            for (int i = 0; i <= N1; i++) file_psi << psi[i][j] << " ";
            file_psi << endl;
        }
    }

    void print() {
        for (int i = 0; i <= N1; i++) {
            for (int j = 0; j <= N2; j++) {
                cout << omega_n[i][j] << "\t";
            }
            cout << endl;
        }
    }


    ~Grid() { cout << "Memory class Grid cleared." << endl; }
};

// Вывод матрицы
void print_A(vector<vector<double>> A) {
    int N = 2;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << A[i][j] << "\t";
        }
        cout << endl;
    }
}

// Вывод вектора 2x1
void print_a(vector<double> a) {
    int N = 2;

    for (int i = 0; i < N; i++) {
        cout << a[i] << "\t";
    }
}

// Обратная матрица
vector<vector<double>> inverse_A(vector<vector<double>> A) {
    int N = 2;
    double d = A[0][0] * A[1][1] - A[1][0] * A[0][1];
    vector<vector<double>> Am1(N, vector<double>(N));

    Am1[0][0] = A[1][1] / d;
    Am1[0][1] = -A[0][1] / d;
    Am1[1][0] = -A[1][0] / d;
    Am1[1][1] = A[0][0] / d;

    return Am1;
}

// Сложение матриц
vector<vector<double>> sum_AB(vector<vector<double>> A, vector<vector<double>> B) {
    int N = 2;
    vector<vector<double>> C(N, vector<double>(N));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }

    return C;
}

// Вычитание матриц
vector<vector<double>> dif_AB(vector<vector<double>> A, vector<vector<double>> B) {
    int N = 2;
    vector<vector<double>> C(N, vector<double>(N));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }

    return C;
}

// Произведение матриц
vector<vector<double>> mult_AB(vector<vector<double>> A, vector<vector<double>> B) {
    int N = 2;
    vector<vector<double>> C(N, vector<double>(N));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = 0;
            for (int k = 0; k < N; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}

// Сложение векторов
vector<double> sum_ab(vector<double> a, vector<double> b) {
    int N = 2;
    vector<double> c(N);

    for (int i = 0; i < N; i++) {
        c[i] = a[i] + b[i];
    }

    return c;
}

// Произведение матрицы на вектор
vector<double> mult_Ab(vector<vector<double>> A, vector<double> b) {
    int N = 2;
    vector<double> c(N);

    for (int i = 0; i < N; i++) {
        c[i] = 0;
        for (int j = 0; j < N; j++) {
            c[i] += A[i][j] * b[j];
        }
    }

    return c;
}


int main(int argc, char* argv[]) {

    setlocale(LC_ALL, "Russian");

    VesselConstriction vc(strtod(argv[5], nullptr), strtod(argv[6], nullptr));
    VesselValve vv(strtod(argv[7], nullptr), strtod(argv[8], nullptr), strtod(argv[9], nullptr));
    Grid grid(strtod(argv[1], nullptr), strtod(argv[2], nullptr), atoi(argv[3]), atoi(argv[4]), strtod(argv[10], nullptr), strtod(argv[11], nullptr), strtod(argv[12], nullptr), vc, vv);

    double delta = strtod(argv[13], nullptr);
    int i, j, it, it_max = atoi(argv[14]);
    double norma_omega, norma_psi, norma_max;
    double h1h1 = grid.get_h1() * grid.get_h1();
    double h2h2 = grid.get_h2() * grid.get_h2();
    double h1h1m2 = h1h1 * 2;
    double h2h2m2 = h2h2 * 2;
    double h1h2m8 = grid.get_h1() * grid.get_h2() * 8;
    double Rem2 = grid.get_Re() * 2;
    double Reh1h1 = grid.get_Re() * h1h1;
    double Reh2h2 = grid.get_Re() * h2h2;
    double Reh1h1m2 = Reh1h1 * 2;
    double Reh2h2m2 = Reh2h2 * 2;
    double Reh1h2m8 = grid.get_Re() * h1h2m8;
    double qm05 = grid.get_q() / 2;
    double tau = grid.get_tau();

    auto start = high_resolution_clock::now();
    it = 0;

    do {

        int N1 = grid.get_N1(), N2 = grid.get_N2();

        // Первый дробный шаг
        // Цикл по горизонтальным линиям сетки
        for (j = 1; j < N2; j++) {
            
            vector<vector<vector<double>>> A(N1 + 1);
            vector<vector<vector<double>>> B(N1);
            vector<vector<vector<double>>> C(N1 + 1);
            vector<vector<double>> F(N1 + 1);
            vector<vector<vector<double>>> alpha(N1 + 1);
            vector<vector<double>> beta(N1 + 1);
            vector<vector<double>> Y(N1 + 1);

            B[0] = { {0, -2 / h1h1}, {0, 0} };
            C[0] = { {1, 0}, {0, 1} };
            F[0] = { 2 / h1h1 * grid.psi_P(j) - grid.dU(j), grid.psi_P(j) };

            for (i = 1; i < N1; i++) {
                A[i] = { {tau / Reh1h1m2 + tau / h1h2m8 * (grid.get_psi_n(i - 1, j + 1) - grid.get_psi_n(i - 1, j - 1)), 0}, {0, tau / h1h1m2} };
                B[i] = { {tau / Reh1h1m2 - tau / h1h2m8 * (grid.get_psi_n(i + 1, j + 1) - grid.get_psi_n(i + 1, j - 1)), 0}, {0, tau / h1h1m2 } };
                C[i] = { {tau / Reh1h1 + 1, tau / Rem2 * (grid.get_c1(i, j) + grid.get_c2(i, j))},
                        {-tau / 2, tau / h1h1 + 1} };
                F[i] = { grid.get_omega_n(i, j) + tau / h1h2m8 * (grid.get_omega_n(i, j + 1) * (grid.get_psi_n(i + 1, j + 1) - grid.get_psi_n(i - 1, j + 1)) - grid.get_omega_n(i, j - 1) * (grid.get_psi_n(i + 1, j - 1) -
                        grid.get_psi_n(i - 1, j - 1))) + tau / Rem2 * ((grid.get_omega_n(i, j + 1) - 2 * grid.get_omega_n(i, j) + grid.get_omega_n(i, j - 1)) / h2h2 + grid.get_c1(i, j) * qm05),
                        grid.get_psi_n(i, j) + tau * (grid.get_psi_n(i, j + 1) - 2 * grid.get_psi_n(i, j) + grid.get_psi_n(i, j - 1)) / h2h2m2 };
            }

            A[N1] = { {0, -2 / h1h1}, {0, 0} };
            C[N1] = { {1, 0}, {0, 1} };
            F[N1] = { 2 / h1h1 * grid.psi_P(j) - grid.dU(j), grid.psi_P(j) };

            // Прямой ход прогонки
            alpha[1] = mult_AB(inverse_A(C[0]), B[0]);
            beta[1] = mult_Ab(inverse_A(C[0]), F[0]);

            for (i = 1; i < N1; i++) {
                alpha[i + 1] = mult_AB(inverse_A(dif_AB(C[i], mult_AB(A[i], alpha[i]))), B[i]);
                beta[i + 1] = mult_Ab(inverse_A(dif_AB(C[i], mult_AB(A[i], alpha[i]))), sum_ab(mult_Ab(A[i], beta[i]), F[i]));
            }

            // Обратный ход прогонки
            Y[N1] = mult_Ab(inverse_A(dif_AB(C[N1], mult_AB(A[N1], alpha[N1]))), sum_ab(mult_Ab(A[N1], beta[N1]), F[N1]));
            for (i = N1 - 1; i >= 0; i--) Y[i] = sum_ab(mult_Ab(alpha[i + 1], Y[i + 1]), beta[i + 1]);

            for (i = 0; i <= N1; i++) { grid.set_omega_np05(i, j, Y[i][0]); grid.set_psi_np05(i, j, Y[i][1]); }
        }

        // Второй дробный шаг
        // Цикл по вертикальным линиям сетки
        for (i = 1; i < N1; i++) {
            vector<vector<vector<double>>> A(N2 + 1);
            vector<vector<vector<double>>> B(N2);
            vector<vector<vector<double>>> C(N2 + 1);
            vector<vector<double>> F(N2 + 1);
            vector<vector<vector<double>>> alpha(N2 + 1);
            vector<vector<double>> beta(N2 + 1);
            vector<vector<double>> Y(N2 + 1);

            B[0] = { {0, -2 / h2h2}, {0, 0} };
            C[0] = { {1, 0}, {0, 1} };
            F[0] = { 0, 0 };

            for (j = 1; j < N2; j++) {
                A[j] = { {tau / Reh2h2m2 - tau / h1h2m8 * (grid.get_psi_np05(i + 1, j - 1) - grid.get_psi_np05(i - 1, j - 1)), 0}, {0, tau / h2h2m2} };
                B[j] = { {tau / Reh2h2m2 + tau / h1h2m8 * (grid.get_psi_np05(i + 1, j + 1) - grid.get_psi_np05(i - 1, j + 1)), 0}, {0, tau / h2h2m2} };
                C[j] = { {tau / Reh2h2 + 1, tau / Rem2 * (grid.get_c1(i, j) + grid.get_c2(i, j))}, {-tau / 2, tau / h2h2 + 1} };
                F[j] = { grid.get_omega_np05(i, j) - tau / h1h2m8 * (grid.get_omega_np05(i + 1, j) * (grid.get_psi_np05(i + 1, j + 1) - grid.get_psi_np05(i + 1, j - 1)) - grid.get_omega_np05(i - 1, j) * (grid.get_psi_np05(i - 1, j + 1) -
                        grid.get_psi_np05(i - 1, j - 1))) + tau / Rem2 * ((grid.get_omega_np05(i + 1, j) - 2 * grid.get_omega_np05(i, j) + grid.get_omega_np05(i - 1, j)) / h1h1 + grid.get_c1(i, j) * qm05),
                        grid.get_psi_np05(i, j) + tau / h1h1m2 * (grid.get_psi_np05(i + 1, j) - 2 * grid.get_psi_np05(i, j) + grid.get_psi_np05(i - 1, j)) };
            }

            A[N2] = { {0, -2 / h2h2}, {0, 0} };
            C[N2] = { {1, 0}, {0, 1} };
            F[N2] = { grid.get_q() / h2h2, qm05 };

            // Прямой ход прогонки
            alpha[1] = mult_AB(inverse_A(C[0]), B[0]);
            beta[1] = mult_Ab(inverse_A(C[0]), F[0]);

            for (j = 1; j < N2; j++) {
                alpha[j + 1] = mult_AB(inverse_A(dif_AB(C[j], mult_AB(A[j], alpha[j]))), B[j]);
                beta[j + 1] = mult_Ab(inverse_A(dif_AB(C[j], mult_AB(A[j], alpha[j]))), sum_ab(mult_Ab(A[j], beta[j]), F[j]));
            }

            // Обратный ход прогонки
            Y[N2] = mult_Ab(inverse_A(dif_AB(C[N2], mult_AB(A[N2], alpha[N2]))), sum_ab(mult_Ab(A[N2], beta[N2]), F[N2]));
            for (j = N2 - 1; j >= 0; j--) Y[j] = sum_ab(mult_Ab(alpha[j + 1], Y[j + 1]), beta[j + 1]);

            for (j = 0; j < N2 + 1; j++) { grid.set_omega_np1(i, j, Y[j][0]); grid.set_psi_np1(i, j, Y[j][1]); }
        }

        for (j = 1; j < N2; j++) {
            grid.set_omega_np1(0, j, 2 / h1h1 * (grid.psi_P(j) - grid.get_psi_np1(1, j)) - grid.dU(j));
            grid.set_omega_np1(N1, j, 2 / h1h1 * (grid.psi_P(j) - grid.get_psi_np1(N1 - 1, j)) - grid.dU(j));
        }

        norma_omega = 0;
        norma_psi = 0;

        for (i = 0; i <= N1; i++) {
            for (j = 0; j <= N2; j++) {
                norma_omega += (grid.get_omega_np1(i, j) - grid.get_omega_n(i, j)) * (grid.get_omega_np1(i, j) - grid.get_omega_n(i, j));
                norma_psi += (grid.get_psi_np1(i, j) - grid.get_psi_n(i, j)) * (grid.get_psi_np1(i, j) - grid.get_psi_n(i, j));
            }
        }

        norma_max = max(sqrt(norma_omega), sqrt(norma_psi));

        for (i = 0; i <= N1; i++) {
            for (j = 0; j <= N2; j++) {
                grid.set_omega_n(i, j, grid.get_omega_np1(i, j)); grid.set_psi_n(i, j, grid.get_psi_np1(i, j));
            }
        }

        it += 1;

    } while (norma_max > delta && it < it_max);

    auto stop = high_resolution_clock::now();
    auto time = duration_cast<seconds>(stop - start);
    grid.saveInFile(it, norma_max, delta, time.count(), vc, vv);

    return 0;
}
