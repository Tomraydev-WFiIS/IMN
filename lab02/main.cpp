#include <cmath>
#include <iostream>
#include <string>

// • Przyjąć następujące wartości parametrów: β = 0.001, N = 500, γ = 0.1, tmax = 100,
// ∆t = 0.1, u0 = 1 (co najmniej jeden osobnik musi być zarażony), TOL = 10−6 , µ <= 20.

#define beta 0.001
#define N 500.0
#define gamma 0.1
#define tmax 100.0
#define delta_t 0.1
#define TOL 1.0e-6
#define mu_max 20

double f(double u)
{
    return (beta * N - gamma) * u - beta * u * u;
}

double df(double u)
{
    double alfa = beta * N - gamma;
    return (alfa - 2 * beta * u);
}

double F1(double U1, double U2, double u_n, double a11, double a12)
{
    double alfa = beta * N - gamma;
    return (U1 - u_n - delta_t * (a11 * (alfa * U1 - beta * U1 * U1) + a12 * (alfa * U2 - beta * U2 * U2)));
}
double F2(double U1, double U2, double u_n, double a21, double a22)
{
    double alfa = beta * N - gamma;
    return (U2 - u_n - delta_t * (a21 * (alfa * U1 - beta * U1 * U1) + a22 * (alfa * U2 - beta * U2 * U2)));
}

double m11(double U1, double a11)
{
    double alfa = beta * N - gamma;
    return (1 - delta_t * a11 * (alfa - 2.0 * beta * U1));
}

double m12(double U2, double a12)
{
    double alfa = beta * N - gamma;
    return (-delta_t * a12 * (alfa - 2.0 * beta * U2));
}

double m21(double U1, double a21)
{
    double alfa = beta * N - gamma;
    return (-delta_t * a21 * (alfa - 2.0 * beta * U1));
}

double m22(double U2, double a22)
{
    double alfa = beta * N - gamma;
    return (1 - delta_t * a22 * (alfa - 2.0 * beta * U2));
}

int main(void)
{
    // 1.1
    // • Rozwiązać równanie (2) stosując metodę trapezów z iteracją Picarda (wzór 10).
    // Sporządzić wykresy u(t) oraz z(t) = N − u(t) i umieścić je na jednym rysunku. (25 pkt)

    // u - zakazeni, z - zdrowi
    int t_n = int(tmax / delta_t);
    double u[t_n];
    double z[t_n];
    double u0 = 1.0;
    u[0] = u0;
    z[0] = N - u0;
    for (int n = 0; n < t_n; n++) {
        double u_mu = u[n]; // punkt startowy
        double alfa = beta * N - gamma;
        for (int mu = 0; mu < mu_max; mu++) {
            u[n+1] = u[n] + (delta_t / 2.0) * ((alfa * u[n] - beta * u[n] * u[n]) + (alfa * u_mu - beta * u_mu * u_mu));
            if (fabs(u[n + 1] - u_mu) < TOL) {
                u_mu = u[n+1];
                break;
            }
            u_mu = u[n+1];
        }
        z[n+1] = N - u[n+1];
    }
    FILE* f11 = fopen("zad1.1.dat", "w");
    for (int n = 0; n < t_n; n++) {
        fprintf(f11, "%g\t%g\t%g\t\n", n * delta_t, u[n], z[n]);
    }
    fclose(f11);

    // 1.2
    // • Rozwiązać równanie (2) stosując metodę trapezów z iteracją Newtona (wzór 13).
    // Sporządzić wykresy u(t) oraz z(t) = N − u(t) i umieścić je na jednym rysnku. (25 pkt)
    FILE* f12 = fopen("zad1.2.dat", "w");
    double t = 0.0;
    double u_n = u0;
    double u_n1 = u0;
    double u_mu = 0.0;
    int mu = 0;
    //maksymalna nr iteracji w petli wewnetrznej
    int max_iter = 0;

    // pu_nkt startowy
    fprintf(f12, "%g %g %g \n", t, u_n1, N - u_n1);
    for (t = delta_t; t < tmax; t += delta_t) {
        mu = 0;
        u_n = u_n1;
        u_mu = 0.;
        while (fabs(u_n1 - u_mu) > TOL && mu < mu_max) {
            u_mu = u_n1;
            u_n1 = u_mu - (u_mu - u_n - (delta_t / 2.0) * (f(u_n) + f(u_mu))) / (1 - (delta_t / 2.0) * df(u_mu));
            if (mu > max_iter) {
                max_iter = mu;
            }
            mu++;
        }
        fprintf(f12, "%g %g %g \n", t, u_n1, N - u_n1);
    }
    fclose(f12);

    // 2
    // Przyjąć następujące wartości parametrów: β = 0.001, N = 500, γ = 0.1, tmax = 100, ∆t = 0.1,
    // u0 = 1 (co najmniej jeden osobnik musi być zarażony), TOL = 10−6, µ ¬ 20.
    // Tak samo jak w pierwszym zadaniu

    // • Rozwiązać równanie (2) stosując metodę niejawną metodę RK, tj. zastosować wzór korektora
    // (17) w którym U1 i U2 wyznaczane są w każdym kroku iteracyjnie według wzorów (20, 21) oraz
    // (27, 28). Sporządzić wykresy u(t) oraz z(t) = N −u(t) i umieścić je na jednym rysnku. (50 pkt)

    FILE* f2 = fopen("zad2.dat", "w");
    double a11 = 0.25;
    double a12 = 0.25 - sqrt(3) / 6;
    double a21 = 0.25 + sqrt(3) / 6;
    double a22 = 0.25;

    double b1 = 0.5;
    double b2 = 0.5;

    double U1 = 0.0; // warnkek startowy
    double U2 = 0.0; // warnkek startowy
    double U1mi = 0.0;
    double U2mi = 0.0;

    t = 0.0;
    u0 = 0.0;
    u_n = u0;
    u_n1 = u0;
    u_mu = 0.0;
    mu = 0;
    fprintf(f2, "%g %g %g \n", t, u_n1, N - u_n1); // punkt startowy


    // u_n1 = u_n + delta_t * (b1*f(U1) * b2*f(U2));

    fclose(f2);
}
