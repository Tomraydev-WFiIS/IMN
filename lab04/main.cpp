#include <stdlib.h>
#include <cmath>
#include <ctime>    
#include <iostream>

#define eps 1.0
#define delta 0.1
#define nx 150
#define ny 100
#define V1 10.0   // top boundary condition
#define V2 0.0    // bottom boundary condition
#define x_max (delta)*(nx)
#define y_max (delta)*(ny)
#define sigma_x (0.1)*(x_max)
#define sigma_y (0.1)*(y_max)
#define TOL 1e-8

double rho[nx+1][ny+1]; // charge density map

double rho1 (int i, int j){
    return exp(
        -pow((i*delta - 0.35*x_max),2)/pow(sigma_x,2)
        -pow((j*delta - 0.5*y_max),2)/pow(sigma_y,2)
    );
}

double rho2 (int i, int j){
    return -exp(
        -pow((i*delta - 0.65*x_max),2)/pow(sigma_x,2)
        -pow((j*delta - 0.5*y_max),2)/pow(sigma_y,2)
    );
}

double sum(double V[][ny + 1] ){
    double S = 0.0;
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            S += delta*delta*(
                0.5*pow((V[i+1][j] - V[i][j])/delta,2) +
                0.5*pow((V[i][j+1] - V[i][j])/delta,2) -
                rho[i][j]*V[i][j]
            );
        }
    }
    return S;
}



void poisson_global(FILE * f_integral, FILE * f_map, FILE * f_err, double omega_g){
    double Vs[nx+1][ny+1] = {0.0}; // old potential
    double Vn[nx+1][ny+1] = {0.0}; // new potential

    // boundary conditions
    for(int i = 0; i <= nx; i++){
        Vs[i][0] = V1;
        Vs[i][ny] = V2;
        Vn[i][0] = V1;
        Vn[i][ny] = V2;
    }
    double S1 = sum(Vs);
    double S = S1;
    double diff;

    std::clock_t start = std::clock();
    double duration;
    int it = 0;
    do
    {
        // Calculate all values except those at the boundary.
        for (int i = 1; i < nx; i++){
            for (int j = 1; j < ny; j++){
                Vn[i][j] = 0.25 * (Vs[i+1][j] + Vs[i-1][j] + Vs[i][j+1] + Vs[i][j-1] + delta*delta/eps*rho[i][j]);
            }
        }

        // Apply Von Neumann boundary condition.
        for (int j = 1; j < ny; j++){
            Vn[0][j] = Vn[1][j];
            Vn[nx][j] = Vn[nx-1][j];
        }

        // Mix the two solutions.
        for (int i = 0; i <= nx; i++){
            for (int j = 0; j <= ny; j++){
                Vs[i][j] = (1.0 - omega_g) * Vs[i][j] + omega_g*Vn[i][j]; // omega_g ∈ (0,1]
            }
        }

        // Calculate the stop condition.
        S = S1;
        S1 = sum(Vs);
        fprintf(f_integral, "%d %lf \n", it, S1);
        diff = (fabs(S1 - S) / S);
    } while ( diff > TOL && ++it);

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    printf("omega_g = %lf\titerations = %d\tduration=%lf s\n", omega_g, it, duration);

    // Save the solution
    for (int i = 0; i <= nx; i++){
        for (int j = 0; j <= ny; j++){
            fprintf(f_map, "%lf %lf %lf \n", delta*i, delta*j, Vs[i][j]); // potential map

            // Calculate the error
            double err = 0;
            if(i > 0 && j > 0 && i < nx && j < ny){
                err = 
                    (Vs[i+1][j] -2*Vs[i][j] + Vs[i-1][j])/(delta*delta) +
                    (Vs[i][j+1] -2*Vs[i][j] + Vs[i][j-1])/(delta*delta) + rho[i][j]/eps;
                fprintf(f_err, "%lf %lf %lf \n", delta*i, delta*j, err);
            }
            else {
                fprintf(f_err, "%lf %lf %lf \n", delta*i, delta*j, 0.0);
            }
        }
        fprintf(f_map, "\n");
        fprintf(f_err, "\n");
    }
}

void poisson_local(FILE * f_integral, double omega_l){
    double V[nx+1][ny+1] = {0.0};

    // boundary conditions
    for(int i = 0; i <= nx; i++){
        V[i][0] = V1;
        V[i][ny] = V2;
    }
    double S1 = sum(V);
    double S = S1;
    double diff;

    std::clock_t start = std::clock();
    double duration;
    int it = 0;
    do
    {
        // Calculate all values except those at the boundary.
        for (int i = 1; i < nx; i++){
            for (int j = 1; j < ny; j++){
                V[i][j] = (1.0 - omega_l)*V[i][j] + (omega_l/4.0)*(V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + delta*delta/eps*rho[i][j]);
            }
        }

        // Apply Von Neumann boundary condition.
        for (int j = 1; j < ny; j++){
            V[0][j] = V[1][j];
            V[nx][j] = V[nx-1][j];
        }

        // Calculate the stop condition.
        S = S1;
        S1 = sum(V);
        fprintf(f_integral, "%d %lf \n", it, S1);
        diff = (fabs(S1 - S) / S);
    } while ( diff > TOL && ++it);

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    printf("omega_l = %lf\titerations = %d\tduration=%lf s\n", omega_l, it, duration);
}

int main(void){
    // setting the charge density map
    FILE * f_charge = fopen("charge.dat", "w");
    for (int i = 0; i <= nx; i++){
        for (int j = 0; j <= ny; j++){
            rho[i][j] = rho1(i,j) + rho2(i,j);
            fprintf(f_charge, "%lf %lf %lf \n", delta*i, delta*j, rho[i][j]);
        }
        fprintf(f_charge, "\n");
    }
    fclose(f_charge);

// 2. Rozwiązać równanie Poissona metodą relaksacji globalnej dla ωG = 0.6, 1.0. W obu przypadkach
// na starcie przyjąć V = 0 w całym obszarze (poza górnym i dolnym brzegiem). Jako warunek
// stopu wykorzystać równania (17) i (18) z parametrem T OL = 10−8.
// Na jedym rysunku umieścić wykresy zmian całki funkcjonalnej S = S(it) dla obu przypadków. (30 pkt)
// Narysować mapę zrelaksowanego potencjału V (x, y) oraz błędu rozwiązania δ = ∇2V (x, y) + ρ(x, y)/ε. (30 pkt)
    // omega = 0.6
    FILE * f_integral = fopen("integral_glob_0.6.dat", "w");
    FILE * f_map = fopen("map_glob_0.6.dat", "w");
    FILE * f_err = fopen("err_glob_0.6.dat", "w");

    poisson_global(f_integral, f_map, f_err, 0.6);
    fclose(f_integral);
    fclose(f_map);
    fclose(f_err);

    f_integral = fopen("integral_glob_1.0.dat", "w");
    f_map = fopen("map_glob_1.0.dat", "w");
    f_err = fopen("err_glob_1.0.dat", "w");

    // omega = 1.0
    poisson_global(f_integral, f_map, f_err, 1.0);
    fclose(f_integral);
    fclose(f_map);
    fclose(f_err);

// 3. Rozwiązać równanie Poissona metodą relaksacji lokalnej dla ωL = 1.0, 1.4, 1.8, 1.9. W każdym
// przypadku na starcie przyjąć V = 0 w całym obszarze (poza górnym i dolnym brzegiem). Jako
// warunek stopu wykorzystać równania (17) i (18) z parametrem T OL = 10−8.
// Na jednym rysunku umieścić wykresy zmian całki funkcjonalnej S = S(it) dla wszystkich rozważanych przypadków. (40 pkt)
    f_integral = fopen("integral_loc_1.0.dat", "w");
    poisson_local(f_integral, 1.0);
    fclose(f_integral);

    f_integral = fopen("integral_loc_1.4.dat", "w");
    poisson_local(f_integral, 1.4);
    fclose(f_integral);

    f_integral = fopen("integral_loc_1.8.dat", "w");
    poisson_local(f_integral, 1.8);
    fclose(f_integral);

    f_integral = fopen("integral_loc_1.9.dat", "w");
    poisson_local(f_integral, 1.9);
    fclose(f_integral);

    return 0;
}