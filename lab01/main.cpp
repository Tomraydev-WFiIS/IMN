#include<cmath>
#include<iostream>
#include<string>

void solve(double dt, double t_max, double lambda, const char * filename, const char * filename_err, double (*fun)(double, double, double));
double euler(double y, double dt, double lambda);
double RK2(double y, double dt, double lambda);
double RK4(double y, double dt, double lambda);


int main(void){
    double lambda = -1;
    double dt_1 = 0.01;
    double dt_2 = 0.1;
    double dt_3 = 1;
    double t_max = 5.0;

    // Analytical solution
    FILE* f = fopen("dt_a.dat", "w");
    for(int i = 0; i <= t_max/0.01; i++){
        fprintf(f, "%g\t%g\n",i*dt_1, exp(lambda * i*dt_1));
    }
    fclose(f);

    // Numeric solution - Euler
    solve(dt_1, t_max, lambda, "dt_1.dat", "dt_1_err.dat", euler);
    solve(dt_2, t_max, lambda, "dt_2.dat", "dt_2_err.dat", euler);
    solve(dt_3, t_max, lambda, "dt_3.dat", "dt_3_err.dat", euler);

    // Numeric solution - solve
    solve(dt_1, t_max, lambda, "dt_1_RK2.dat", "dt_1_RK2_err.dat", RK2);
    solve(dt_2, t_max, lambda, "dt_2_RK2.dat", "dt_2_RK2_err.dat", RK2);
    solve(dt_3, t_max, lambda, "dt_3_RK2.dat", "dt_3_RK2_err.dat", RK2);

    // Numeric solution - RK4
    solve(dt_1, t_max, lambda, "dt_1_RK4.dat", "dt_1_RK4_err.dat", RK4);
    solve(dt_2, t_max, lambda, "dt_2_RK4.dat", "dt_2_RK4_err.dat", RK4);
    solve(dt_3, t_max, lambda, "dt_3_RK4.dat", "dt_3_RK4_err.dat", RK4);


}

void solve(double dt,
            double t_max,
            double lambda,
            const char * filename,
            const char * filename_err,
            double (*fun)(double, double, double)
            ){
    FILE* f = fopen(filename, "w");
    FILE* f_err = fopen(filename_err, "w");
    int n = t_max / dt;
    double y[n+1];
    y[0] = 1.0;
    fprintf(f, "0\t%g\n", y[0]);
    fprintf(f_err, "0\t0\n");
    for(int i = 1; i < n+1; i++){
        y[i] = fun(y[i-1], dt, lambda);
        fprintf(f, "%g\t%g\n", i*dt, y[i]);
        fprintf(f_err, "%g\t%g\n",i*dt, y[i] - exp(lambda * i*dt));
    }
    fclose(f);
    fclose(f_err);
}


double euler(double y, double dt, double lambda){
    return y + dt * lambda * y;
}

double RK2(double y, double dt, double lambda){
    double k1 = lambda*y;
    double k2 = lambda*(y + dt*k1);
    return y + (dt / 2.0) * (k1 + k2);
}

double RK4(double y, double dt, double lambda){
    double k1 = lambda*y;
    double k2 = lambda*(y + dt/2.0*k1);
    double k3 = lambda*(y + dt/2.0*k2);
    double k4 = lambda*(y + dt*k3);
    return y + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4);
}