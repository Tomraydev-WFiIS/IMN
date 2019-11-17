#include <cstdio>
#include <cstdlib>
#include <cmath>

//predefiniowane stale
#define delta 0.2
#define dx 0.2
#define dy 0.2
#define nx 128
#define ny 128

//stale wyrazenia
const double xmax = delta*nx;
const double ymax = delta*ny;
const double TOL = pow(10,-8);
//warunek stopu
double stop(double V[][ny+1], int k) {
	double s = 0.;
	for(int i = 0; i <= nx-k; i += k){
		for(int j = 0; j <= ny-k; j += k){
			s += 0.5*pow(k*delta, 2)*(pow((V[i+k][j] - V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i][j+k])/(2*k*delta), 2) +
					pow((V[i][j+k] - V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i+k][j])/(2*k*delta), 2));
		}
	}
	return s;
}


//relaksacja wielosiatkowa
void poisson(FILE * calka, FILE * mapy){
    //stopy
    double s = 0.;
    double s1 = 0.;
    int iter = 0;
    int k = 16;

    double V[nx+1][ny+1];

    //zerowanie tablicy
    for(int i = 0 ; i <= nx; i++){
        for(int j = 0; j <= ny; j++){
            V[i][j] = 0.;
        }
    }
    //warunki brzegowe
	for(int j = 0; j <= ny; j++){
		V[0][j] = sin(M_PI*dy*j/ymax);
		V[nx][j] = sin(M_PI*dy*j/ymax);
	}

	for(int i = 0; i <= nx; i++){
		V[i][ny] = (-1)*sin(2*M_PI*dx*i/xmax);
		V[i][0] = sin(2*M_PI*dx*i/xmax);
	}
   

    while(k > 0){
        s1 = stop(V, k);

        do{
			for(int i = k; i <= nx-k; i += k){
		        for(int j = k; j <= ny-k; j += k){
			        V[i][j] = 0.25*(V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);
                }
            }

            s = s1;
			s1 = stop(V, k);
            fprintf(calka, "%d %d %f \n",k, iter, s1 );
            iter++;
		}while(fabs((s1-s)/s) > TOL);

        fprintf(calka, "\n \n" );

        for(int i = 0; i <= nx; i += k){
            for(int j = 0; j <= ny; j += k){
                fprintf(mapy, "%f %f %f \n", dx*i, dy*j, V[i][j]);
            }
            fprintf(mapy, "\n");
        }
        fprintf(mapy, "\n");

	    for(int i = 0; i <= nx-k; i += k){
		    for(int j = 0; j <= ny-k; j += k){
                V[i + k/2][j + k/2] = 0.25*(V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);
                if(i!=nx-k)
                    V[i + k][j + k/2] = 0.5*(V[i+k][j] + V[i+k][j+k]);
                if(j!=ny-k)
                    V[i + k/2][j + k] = 0.5*(V[i][j+k] + V[i+k][j+k]);
                if(j!=0)
                    V[i + k/2][j] = 0.5*(V[i][j] + V[i+k][j]);
                if(i!=0)
                    V[i][j + k/2] = 0.5*(V[i][j] + V[i][j+k]);
		    }
	    }
	    k = k/2;
	}
}

int main(){
    FILE * calka = fopen("calka.dat", "w");
    FILE * mapy = fopen("mapy.dat", "w");
    poisson( calka, mapy);
    fclose(calka);
    fclose(mapy);

    return 0;
}
