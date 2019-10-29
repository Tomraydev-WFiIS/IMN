#include<cmath>
#include<iostream>
#include<string>

double f(int t, int u, double beta, int N, double gamma){
    return (beta * (double)N - gamma)*(double)u - beta*(double)u*(double)u;
}


int main(void){
    // • Przyjąć następujące wartości parametrów: β = 0.001, N = 500, γ = 0.1, tmax = 100,
    // ∆t = 0.1, u0 = 1 (co najmniej jeden osobnik musi być zarażony), TOL = 10−6 , µ <= 20.
    double beta = 0.001;
    int N = 500;
    double gamma = 0.1;
    int tmax = 100;
    double delta_t = 0.1;
    double u0 = 1.0;
    double TOL = 1.0e-6;


    // 1.
    // • Rozwiązać równanie (2) stosując metodę trapezów z iteracją Picarda (wzór 10).
    // Sporządzić wykresy u(t) oraz z(t) = N − u(t) i umieścić je na jednym rysunku. (25 pkt)


    // u - zakazeni, z - zdrowi
    int t_n = int(tmax/delta_t);
    double u[t_n];
    int z[t_n];
    u[0] = u0;
    z[0] = N - u0;
    for(int n = 0; n < t_n; n++){
        double u_mu = u[n]; // punkt startowy
        double alfa = beta * (double)N - gamma;
        for(int mu = 0; mu < 20; mu++){
            u[n+1] = u[n] + (delta_t/2.0) *
                (alfa * u[n] - beta*u[n]*u[n]) +
                (alfa * u_mu - beta*u_mu*u_mu);
            if (abs(u[n+1] - u_mu) < TOL){
                u[n+1] = u_mu;
                break;
            }
            u[n+1] = u_mu;
        }
        printf("%.4f\n", u[n+1]);
    }


    // • Rozwiązać równanie (2) stosując metodę trapezów z iteracją Newtona (wzór 13).
    // Sporządzić wykresy u(t) oraz z(t) = N − u(t) i umieścić je na jednym rysunku. (25 pkt)


    // 2.
    // Przyjąć następujące wartości parametrów: β = 0.001, N = 500, γ = 0.1, tmax = 100, ∆t = 0.1,
    // u0 = 1 (co najmniej jeden osobnik musi być zarażony), TOL = 10−6, µ ¬ 20.

    // • Rozwiązać równanie (2) stosując metodę niejawną metodę RK, tj. zastosować wzór korektora
    // (17) w którym U1 i U2 wyznaczane są w każdym kroku iteracyjnie według wzorów (20, 21) oraz
    // (27, 28). Sporządzić wykresy u(t) oraz z(t) = N −u(t) i umieścić je na jednym rysunku. (50 pkt)


    }



// u[n+1] = u[n] + (delta_t/2.0) * (f(delta_t*n,u[n], beta, N, gamma) + f(delta_t*(n+1),u[n+1], beta, N, gamma));