#include <stdlib.h>
#include <cmath>
#include <ctime>    
#include <iostream>
#include <fstream>

#define delta  0.01 // delta x, delta y
#define nx 400
#define ny 90
#define i1 200
#define i2 210
#define j1 50
#define sigma (10.0*(delta))
#define xa 0.45
#define ya 0.45
#define ITMAX 20


// #define y_j1 ((j1)*(delta))
// #define y_ny ((ny)*(delta))
// #define x_i1 ((i1)*(delta))
// #define x_nx ((nx)*(delta))

void set_velocity_field(double vx[][ny+1], double vy[][ny+1], double psi[][ny+1]){
    //wyznacz mape predkosci
    FILE * f_vx = fopen("vx.dat", "w");
    FILE * f_vy = fopen("vy.dat", "w");
    for (int i = 1; i < nx; i++){
        for (int j = 1; j < ny; j++){
            vx[i][j] = (psi[i][j+1] - psi[i][j-1])/ (2.0*delta);
            vy[i][j] = -(psi[i+1][j] - psi[i-1][j])/ (2.0*delta);
            fprintf(f_vx, "%d %d %g\n", i, j, vx[i][j]);
            fprintf(f_vy, "%d %d %g\n", i, j, vy[i][j]);
        }
        fprintf(f_vx, "\n");
        fprintf(f_vy, "\n");
    }
    fclose(f_vx);
    fclose(f_vy);

    // zastawka
    for (int i = i1; i < i2+1; i++){
        for (int j = 0; j < j1+1; j++){
            vx[i][j] = (psi[i][j+1] - psi[i][j-1])/ (2.0*delta);
            vx[i][j] = -(psi[i][j+1] - psi[i][j-1])/ (2.0*delta);
        }
    }
    // dolny i górny brzeg
    for (int i = 1; i < nx; i++){
        vx[i][0] = vy[i][ny] = 0.0;
    }
    // lewy i prawy brzeg
    for (int j = 0; j < ny+1; j++){
        vx[0][j] = vx[1][j];
        vy[nx][j] = vy[nx-1][j];
    }
}

bool is_inside_valve(int i, int j){
    if (i >= i1 && i <=i2 && j >= 0 && j <= j1){
        return true;
    }else{
        return false;
    }
}

double eq9_wb(int i, int j, double u0[][ny+1], double u1[][ny+1], double vx[][ny+1], double vy[][ny+1], double delta_t){
    double D = 0;
    double l1 = (1.0/(1+ ((2.0*D*delta)/pow(delta,2)))) * (u0[i][j] - delta_t/2.0*vx[i][j] *((u1[i+1][j] - u1[nx][j])/(2.0*delta) + (u0[i+1][j] - u0[nx][j])/(2.0*delta)));
    double l2 = delta_t/2.0 *vy[i][j] * ((u0[i][j+1] - u0[i][j-1]) /(2.0*delta) + (u1[i][j+1] - u1[i][j-1]) /(2.0*delta));
    double l3 = delta_t/2.0*D * ((u0[i+1][j] + u0[nx][j] + u0[i][j+1] + u0[i][j-1] - 4.0*u0[i][j])/pow(delta,2));
    double l4 = (u1[i+1][j] + u1[nx][j] + u1[i][j+1] + u1[i][j-1])/(2.0 * delta);
    return l1-l2+l3+l4;
}

double eq9(int i, int j, double u0[][ny+1], double u1[][ny+1], double vx[][ny+1], double vy[][ny+1], double delta_t){
    //TODO
}

void ad(double u0[][ny+1], double u1[][ny+1], double vx[][ny+1], double vy[][ny+1], double delta_t){
    // algorytm adwekcji dyfuzji
    for (int IT = 0; IT <= ITMAX; IT++){
        // start iteracja Picarda
        for (int i = 0; i < nx+1; i++){
            for (int j = 0; j < ny+1; j++){
                u1[i][j] = u0[i][j];
            }
        }
        for (int k = 1; k <= 20; k++){
            for (int i = 0; i <= nx; i++){
                for (int j = 0; j < ny; j++){
                    if(is_inside_valve(i,j)) {
                        continue;
                    }else if(i==0 || i == nx){
                        u1[i][j] = eq9_wb(i,j, u0, u1, vx, vy, delta_t);
                    }else {
                        u1[i][j] = eq9(i,j, u0, u1, vx, vy, delta_t);
                    }
                }
            }
        }
    } // IT
    for (int i = 0; i < nx+1; i++){
        for (int j = 0; j < ny+1; j++){
            u0[i][j] = u1[i][j];
        }
    }
}
int main(void){
    // initialize map
    double psi[nx+1][ny+1] = {0};
    double vx[nx+1][ny+1] = {0};
    double vy[nx+1][ny+1] = {0};
    double u0[nx+1][ny+1] = {0}; // gestosć
    double u1[nx+1][ny+1] = {0}; // gestosć

    // wczytaj mape psi
    std::ifstream f_psi;
    f_psi.open("psi.dat");
    int x;
    int y;
    double val;
    while (f_psi >> x){
        f_psi >> y;
        f_psi >> val;
        psi[x][y] = val;
    }
    set_velocity_field(vx, vy, psi);
    
    // maksymalna predkosc
    double v_max = 0.0;
    double v = 0.0;
    for (int i = 0; i < nx+1; i++){
        for (int j = 0; j < ny+1; j++){
            v = sqrt(pow(vx[i][j],2)+pow(vy[i][j],2));
            if(v>v_max) {v_max=v;}
        }
    }
    printf("v_max = %g\n", v_max);

    // warunek początkowy
    for (int i = 0; i < nx+1; i++){
        for (int j = 0; j < ny+1; j++){
            u0[i][j] = (1.0/2*M_PI*sigma*sigma)
            *exp(- (x-xa)*(x-xa) + (y-ya)*(y-ya) / (2*sigma*sigma));
        }
    }

    double delta_t = delta/v_max; // krok czasowy
    ad(u0, u1, vx, vy, delta_t);
}