#include <stdio.h>
#include <math.h>
#include "mgmres.h"
#include "mgmres.c"

#define delta 0.1
#define N ((nx+1)*(ny+1))

double rho1 (int i, int j, int nx, int ny){
    double y_max = delta*ny;
    double x_max = delta*nx;
    double sigma = x_max/10.0;
    double x = i*delta;
    double y = j*delta;

    return exp(
        -pow((i*delta - 0.35*x_max),2)/pow(sigma,2)
        -pow((j*delta - 0.5*y_max),2)/pow(sigma,2)
    );
}

double rho2 (int i, int j, int nx, int ny){
    double y_max = delta*ny;
    double x_max = delta*nx;
    double sigma = x_max/10.0;
    double x = i*delta;
    double y = j*delta;

    return -exp(
        -pow((i*delta - 0.65*x_max),2)/pow(sigma,2)
        -pow((j*delta - 0.5*y_max),2)/pow(sigma,2)
    );
}


int j(int l, int nx){
    return floor(l/(nx+1.0));
}

int i(int l, int nx){
    return  l - j(l,nx)*(nx+1.0);
}


double el(int l, int e1, int e2, int nx){
    if (i(l, nx) <= (nx/2)){
        return e1;
    }
    else {
        return e2;
    }
}

int fill_A_matrix(FILE * f_matrix, FILE * f_vector,int nx, int ny,int V1,int V2,int V3,int V4,int e1, int e2,int *ia, int *ja, double *b, double *a){
    // WB Dirichleta zgodnie z algorytmem przedstawionym w sekcji (3).
    // W celu sprawdzenia poprawności wypełnienia macierzy A i wektora b należy zapisać je do pliku, podając niezerowe elementy l, il , jl , a[l] (dla macierzy) oraz l, il , jl , b[l] (dla wektora). 20 pkt.
    int k = -1;
    int brzeg; // wskaźnik położenia : 0 - środek obszaru ; 1 - brzeg
    double vb; // potencjal na brzegu
    for (int l = 0; l < N; l++){ //DO
        brzeg = 0;
        vb = 0;

        if(i(l, nx) == 0){ //lewy brzeg
            brzeg = 1;
            vb = V1;
        }

        if(j(l, nx) == ny){ //górny brzeg
            brzeg = 1;
            vb = V2;
        }

        if(i(l, nx) == nx){ //prawy brzeg
            brzeg = 1;
            vb = V3;
        }

        if(j(l, nx) == 0){ //dolny brzeg
            brzeg = 1;
            vb = V4;
        }

        // wypełniamy od razu wektor wyrazów wolnych
        b[l] = -rho1(i(l,nx), j(l,nx),nx,ny);

        if(brzeg == 1){
            b[l] = vb; // wymuszamy wartość potencjału na brzegu
        }

        // wypełniamy elementy macierzy A
        ia[l] = -1; // wskaźnik dla pierwszego el . w wierszu

        // lewa skrajna przekatna
        if(l-nx-1 >= 0 && brzeg==0){
            k++;
            if(ia[l] < 0){
                ia[l] = k;
            }
            a[k] = el(l,e1,e2,nx) / (delta*delta); // a_l_l-nx-1
            ja[k] = l-nx-1;
        }

        // poddiagonala
        if(l-1 >= 0 && brzeg==0){
            k++;
            if(ia[l] < 0){
                ia[l] = k;
            }
            a[k] = el(l,e1,e2,nx) / (delta*delta); //a_l_l-1
            ja[k] = l-1;
        }

        // diagonala
        k++;
        if(ia[l]<0){
            ia[l] = k;
        }
        if(brzeg==0){
            a[k] = -(2*el(l,e1,e2,nx) + el(l+1,e1,e2,nx) + el(l+nx+1,e1,e2,nx)) / (delta*delta);
        }
        else{
            a[k]=1.0;
        }
        ja[k] = l;

        //naddiagonala
        if(l<N && brzeg ==0){
            k++;
            a[k] = (el(l+1,e1,e2,nx))/(delta*delta); // a_l_l+1
            ja[k] = l+1;
        }

        // prawa skrajna przekątna
        if(l < N-nx-1 && brzeg ==0){
            k++;
            a[k] = el(l+nx+1,e1,e2,nx)/(delta*delta);// a_l_l+nx+1
            ja[k]=l+nx+1;
        }
        fprintf(f_vector, "%d\t%d\t%d\t%lf\n",l, i(l,nx) , j(l,nx), b[l]);
    } //END DO
    int nz_num = k+1;
    ia[N] = nz_num; // ilosc niezerowych elementow (1 element ma indeks 0)
    for(int k = 0; k < 5*N; k++){
        fprintf(f_matrix, "%d\t%lf\n",k, a[k]);
    }
    return nz_num;
}

void solve(FILE * f_map, FILE * f_matrix, FILE * f_vector,int nx, int ny,int V1,int V2,int V3,int V4,int e1, int e2,int *ia, int *ja, double *a){
    double b[N];
    double V[N];
    // 3. Wyznaczamy niezerowe elementy macierzy A oraz wektora wyrazów wolnych b uwzględniając
    int nz_num = fill_A_matrix(f_matrix, f_vector, nx,ny, V1, V2, V3, V4, e1, e2, ia, ja, b, a);
    // 4. Rozwiązujemy układ równań z macierzą rzadką stosując procedurę (dołączyć pliki: mgmres.c i) mgmres.h)
    int itr_max = 500;
    int mr = 500;
    double tol_abs = 1e-8;
    double tol_rel = 1e-8;
    pmgmres_ilu_cr(N,nz_num,ia,ja,a,V,b,itr_max, mr,tol_abs,tol_rel);
    // 5. Sporządzić mapy potencjału:
    double limit = 0.0;
    for(int l = 0; l < N; l++){
        if(delta*j(l,nx) > limit){
            fprintf(f_map,"\n");
        }
        fprintf(f_map,"%lf\t%lf\t%lf\n", delta*j(l, nx), delta*i(l, ny), V[l]);
        limit = delta*j(l, nx);
    }

}

int main(void){
    // 1. Ustalamy początkowe wartości parametrów: ∆ = 0.1, nx = ny = 4, ε1 = ε2 = 1, V1 = V3 = 10,
    // V2 = V4 = −10, ρ1(x, y) = ρx,y = 0, N = (nx + 1)(ny + 1)
    FILE * f_matrix = fopen("matrix.dat", "w");
    FILE * f_vector = fopen("vector.dat", "w");
    FILE * f_map4 = fopen("map4.dat", "w");
    FILE * f_map50 = fopen("map50.dat", "w");
    FILE * f_map100 = fopen("map100.dat", "w");
    FILE * f_map200 = fopen("map200.dat", "w");
    int nx = 4;
    int ny = 4;
    int e1 = 1;
    int e2 = 1;
    int V1 = 10;
    int V2 = -10;
    int V3 = 10;
    int V4 = -10;
    double rho0 = 0.0;

    double a[5*N]; //niezerowe wartości el. macierzowych
    int ja[5*N]; // (przechowuje informacje o numerach kolumn
    int ia[N+1]; //wskaźniki do elementów rozpoczynających dany wiersz
    for(int i = 0;i < N+1; i++){
        ia[i] = -1;
    }
    solve(f_map4, f_matrix, f_vector, nx,ny, V1, V2, V3, V4, e1, e2, ia, ja, a);
    fclose(f_map4);
    //a) nx = ny = 50,
    nx = 50;
    ny = 50;
    solve(f_map50, f_matrix, f_vector, nx,ny, V1, V2, V3, V4, e1, e2, ia, ja, a);
    fclose(f_map50);
    //b) nx = ny = 100,
    nx = 100;
    ny = 100;
    solve(f_map100, f_matrix, f_vector, nx,ny, V1, V2, V3, V4, e1, e2, ia, ja, a);
    fclose(f_map100);

    //c) nx = ny = 200
    nx = 200;
    ny = 200;
    solve(f_map200, f_matrix, f_vector, nx,ny, V1, V2, V3, V4, e1, e2, ia, ja, a);
    fclose(f_map200);

    // 6. Znaleźć rozkłady potencjału dla 3 przypadków:
    //  a) ε1 = 1 i ε2 = 1, b) ε1 = 1 i ε2 = 2 oraz c) ε1 = 1 i ε2 = 10
}
