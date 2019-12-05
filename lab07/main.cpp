#include <stdlib.h>
#include <cmath>
#include <ctime>    
#include <iostream>

#define delta  0.01 // delta x, delta y
#define rho  1.0   // gestosc
#define mu  1.0    // lepkosc
#define nx  200
#define ny  90
#define i1  50
#define j1  55
#define IT_MAX  20000

#define y_j1 ((j1)*(delta))
#define y_ny ((ny)*(delta))
#define x_i1 ((i1)*(delta))
#define x_nx ((nx)*(delta))

// zeta - funkcja wirowosci
// psi - funkcja strumienia

void set_psi_boundary(double psi[][ny+1], double Q_we){
    double Q_wy = Q_we * (pow(y_ny,3) - pow(y_j1,3) - 3*y_j1*y_ny*y_ny + 3*y_j1*y_j1*y_ny) / (pow(y_ny,3));
    double x,y;
    // A - wejscie
    for (int j = j1; j <= ny; j++){
        y = delta*j;
        psi[0][j] = Q_we/(2*mu) * (pow(y,3)/3.0 - pow(y,2)/2.0 * (y_j1+y_ny) + y*y_j1*y_ny);
    }

    // C - wyjscie
    for (int j = 0; j <= ny; j++){
        y = delta*j;
        psi[nx][j] = Q_wy/(2*mu) * (pow(y,3)/3.0 - pow(y,2)/2.0*y_ny)
            + (Q_we * pow(y_j1,2)*(-y_j1 + 3.0*y_ny))/(12*mu);
    }
    
    // B - góra
    for (int i = 1; i < nx; i++){
        psi[i][ny] = psi[0][ny];
    }
    
    // D - dół
    for (int i = i1; i < nx; i++){
        psi[i][0] = psi[0][j1];
    }
    
    // E - prawy bok prostokąta
    for (int j = 1; j <= j1; j++){
        psi[i1][j] = psi[0][j1];
    }
    
    // F - góra prostokąta
    for (int i = 1; i <= i1; i++){
        psi[i][j1] = psi[0][j1];
    }
}

void set_zeta_boundary(double zeta[][ny+1], double psi[][ny+1], double Q_we){
    double Q_wy = Q_we * (pow(y_ny,3) - pow(y_j1,3) - 3*y_j1*y_ny*y_ny + 3*y_j1*y_j1*y_ny) / (pow(y_ny,3));
    double x,y;
    // A - wejscie
    for (int j = j1; j <= ny; j++){
        y = delta*j;
        zeta[0][j] = Q_we/(2.0*mu)*(2.0*y-y_j1-y_ny);
    }

    // C - wyjscie
    for (int j = 0; j <= ny; j++){
        y = delta*j;
        zeta[nx][j] = Q_wy/(2.0*mu)*(2.0*y-y_ny);
    }
    
    // B - góra
    for (int i = 1; i < nx; i++){
        zeta[i][ny] = 2.0/(delta*delta)*(psi[i][ny-1] - psi[i][ny]);
    }
    
    // D - dół
    for (int i = i1+1; i < nx; i++){
        zeta[i][0] = 2.0/(delta*delta)*(psi[i][1] - psi[i][0]);
    }
    
    // E - prawy bok prostokąta
    for (int j = 1; j < j1; j++){
        zeta[i1][j] = 2.0/(delta*delta)*(psi[i1+1][j] - psi[i1][j]);
    }
    
    // F - góra prostokąta
    for (int i = 1; i <= i1; i++){
        zeta[i][j1] = 2.0/(delta*delta)*(psi[i][j1+1] - psi[i][j1]);
    }

    // E/F - wierzcholek
    zeta[i1][j1] = 0.5*(zeta[i1-1][j1] + zeta[i1][j1-1]);
}

void relax(double psi[][ny+1], double zeta[][ny+1], double Q_we){
    double omega;
    double gamma;
    double x,y;
    int j2 = j1+2;
    for (int IT = 1; IT <= IT_MAX; IT++){
        if (IT < 2000){
            omega = 0.0;
        }else{
            omega = 1.0;
        }

        for (int i = 1; i < nx; i++){
            x = delta*i;
            for (int j = 1; j < ny; j++){
                y = delta*j;
                if( (i <= i1 && j > j1) || (i > i1) ){
                    psi[i][j] = 0.25*(psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - delta*delta*zeta[i][j]);
                    zeta[i][j] = 0.25*(zeta[i+1][j] + zeta[i-1][j] + zeta[i][j+1] + zeta[i][j-1]);
                    zeta[i][j] -= omega*rho/(16*mu)*(
                        (psi[i][j+1] - psi[i][j-1])*(zeta[i+1][j]-zeta[i-1][j])
                        -(psi[i+1][j] - psi[i-1][j])*(zeta[i][j+1] - zeta[i][j-1])
                    );
                }
            } // j
        } // i
        // modyfikacja WB zeta
        set_zeta_boundary(zeta, psi, Q_we);
        // kontrola bledu gamma
        gamma = 0;
        for (int i = 1; i < nx; i++){
            gamma += psi[i+1][j2] + psi[i-1][j2]  + psi[i][j2+1]  + psi[i][j2-1] - 4*psi[i][j2] - delta*delta*zeta[i][j2];
        }
        // printf("IT = %d\tΓ = %g\n", IT, gamma); // debug
    } // IT
    
}

void save_file(double psi[][ny+1], double zeta[][ny+1],FILE *f){
    double x, y, u, v;
    for (int i = 0; i < nx+1; i++){
        x = delta*i;
        for (int j = 0; j < ny+1; j++){
            y = delta*j;
            if( (i > 0 && i <= i1 && j > j1 && j < ny) || (i > i1 && i < nx && j > 0 && j < ny) ){
                u = (psi[i][j+1] - psi[i][j-1]) / (2*delta);
                v = (psi[i+1][j] - psi[i-1][j]) / (2*delta);
            }else {
                u = v = 0.0;
            }
            fprintf(f, "%lf %lf %g %g %g %g\n", x, y, psi[i][j], zeta[i][j],u, v);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

int main(void){
    // Ad. 4, Q = -1000
    double psi[nx+1][ny+1] = {0};
    double zeta[nx+1][ny+1] = {0};
    double Q = -1000.0;
    set_psi_boundary(psi,Q);
    set_zeta_boundary(zeta, psi, Q);
    relax(psi, zeta, Q);
    FILE * f = fopen("Q_-1000.txt", "w");
    save_file(psi, zeta, f);
    fclose(f);

    // Ad. 5, Q = -4000
    psi[nx+1][ny+1] = {0};
    zeta[nx+1][ny+1] = {0};
    Q = -4000.0;
    set_psi_boundary(psi,Q);
    set_zeta_boundary(zeta, psi, Q);
    relax(psi, zeta, Q);
    f = fopen("Q_-4000.txt", "w");
    save_file(psi, zeta, f);
    fclose(f);


    // Ad. 6, Q = 4000
    psi[nx+1][ny+1] = {0};
    zeta[nx+1][ny+1] = {0};
    Q = 4000.0;
    set_psi_boundary(psi,Q);
    set_zeta_boundary(zeta, psi, Q);
    relax(psi, zeta, Q);
    f = fopen("Q_+4000.txt", "w");
    save_file(psi, zeta, f);
    fclose(f);
}