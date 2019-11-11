#include <cstdio>
#include <cstdlib>
#include <cmath>

//predefiniowane stale
#define eps 1.
#define delta 0.1
#define dx 0.1
#define dy 0.1
#define nx 150
#define ny 100
#define V1 10
#define V2 0

//stale wyrazenia
const double xmax = delta*nx;
const double ymax = delta*ny;
const double sigX = 0.1*xmax;
const double sigY = 0.1*ymax;
const double TOL = pow(10,-8);

//tablica gestosci globalna
double ro[nx + 1][ny + 1];

//ustawienie tablicy gestosci
void setRo(){
    for(int i = 0; i <= nx; i++){
        for(int j = 0; j <= ny; j++){
			ro[i][j] =  exp( -pow(dx*i - 0.35 * xmax, 2) / (sigX*sigX) - pow(dy*j - 0.5 * ymax, 2) / (sigY*sigY) )
                         + ( -exp( -pow(dx*i - 0.65 * xmax, 2) / (sigX*sigX) - pow(dy*j - 0.5 * ymax, 2) / (sigY*sigY) ) );
		}
	}
}

//warunek stopu petli w algorytmie
double stop(double V[][ny + 1] ){
    double s = 0.;
    	for(int i = 0; i < nx; i++) {
			for(int j = 0; j < ny; j++) {
				s += (delta*delta) * ( 0.5*pow((V[i+1][j] - V[i][j])/delta, 2) +  0.5*pow((V[i][j+1] - V[i][j])/delta, 2) - ro[i][j]*V[i][j]);
			}
		}
    // printf("%g \n",s);

    return s;
}



//relaksacja globalna
void relGlobal(double wg, FILE * calka, FILE * mapa, FILE * blad){

//stopy
    double s = 0.;
    double s1 = 0.;
    int iter = 0;
//tablice
    double Vn[nx+1][ny+1];
    double Vs[nx+1][ny+1];


    //zerowanie Vnowa
    for(int i = 0 ; i <= nx; i++){
        for(int j = 0; j <= ny; j++){
            Vn[i][j] = 0.;
        }
    }
    //warunki brzegowe Vnowa
    for(int i = 0; i <= nx; i++){
        Vn[i][0] = V1;
        Vn[i][ny] = V2;
    }

    //tablica kopia Vstara
    for(int i = 0 ; i <= nx; i++){
        for(int j = 0; j <= ny; j++){
            Vs[i][j] = Vn[i][j];
        }
    }

//algorytm
    s1 = stop(Vn);
    printf("%f", s1);

    do{
        iter++;
        //1 etap
        for(int i = 1; i < nx; i++){
            for(int j = 1; j < ny; j++){
                Vn[i][j] = 0.25 * ( Vs[i+1][j] + Vs[i-1][j] + Vs[i][j+1] + Vs[i][j-1]
                    +  delta*delta/eps * ro[i][j] );
                }
        }
        //2 etap
            for(int j = 1; j < ny; j++){
                Vn[0][j] = Vn[1][j];
                Vn[nx][j] = Vn[nx - 1][j];
            }

        //3 etap
        for(int i = 0 ; i <= nx ; i++){
            for(int j = 0; j <= ny; j++){
                Vs[i][j] = (1-wg)*Vs[i][j] + wg*Vn[i][j];
            }
        }
        s = s1;
        s1 = stop(Vs); //sprawdzic czy to
        //zapis do calki funkcjonalnej
        fprintf(calka, "%d %f \n", iter, s1);
    }while(fabs( (s1 - s)/s ) > TOL);


    for(int i = 0; i <= nx; i++){
        for(int j = 0; j <= ny; j++){
            //zapis do potencjalu
            fprintf(mapa, "%f %f %f \n", dx*i, dy*j, Vs[i][j]);

            //zapis do bledu
            if(i >0 && j > 0 && i < nx && j < ny)
                fprintf(blad, "%f %f %f \n", dx*i, dy*j, (Vs[i+1][j] -2*Vs[i][j] + Vs[i-1][j])/(delta*delta) +
                    (Vs[i][j+1] -2*Vs[i][j] + Vs[i][j-1])/(delta*delta) + ro[i][j]/eps );
             else
                fprintf(blad, "%f %f %f \n", dx*i, dy*j, 0.0);
        }
        fprintf(mapa, "\n");
        fprintf(blad, "\n");

    }

    printf("%d \n", iter);

}

void relLocal(double wl,  FILE * calka){

//stopy
    double s = 0.;
    double s1 = 0.;
    int iter = 0;
//tablice
    double V[nx+1][ny+1];


    //zerowanie Vnowa
    for(int i = 0 ; i <= nx; i++){
        for(int j = 0; j <= ny; j++){
            V[i][j] = 0.;
        }
    }
    //warunki brzegowe Vnowa
    for(int i = 0; i <= nx; i++){
        V[i][0] = V1;
        V[i][ny] = V2;
    }

//algorytm
    s1 = stop(V);
    do{
        iter++;

        //1 etap
        for(int i = 1; i < nx; i++){
            for(int j = 1; j < ny; j++){
                V[i][j] = (1 - wl)* V[i][j] +(wl/4) * ( V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1]
                    +  delta*delta/eps * ro[i][j] );
                }
        }
        //2 etap
            for(int j = 1; j < ny; j++){
                V[0][j] = V[1][j];
                V[nx][j] = V[nx - 1][j];
            }

        s = s1;
        s1 = stop(V);
        //zapis do calki funkcjonalnej
        fprintf(calka, "%d %f \n", iter, s1);
    }while(fabs( (s1 - s)/s ) > TOL);


    printf("%d \n", iter );
}


int main(){

//ustawienie tablicy gestosci
    setRo();
//zad 1
    double wg;

    wg = 1.0;
    FILE * calka = fopen("calka1.0_glob.dat", "w");
    FILE * mapa = fopen("mapa1.0_glob.dat", "w");
    FILE * blad = fopen("blad1.0_glob.dat", "w");
    printf("Relaksacja globalna 1.0, iteracje ");
    relGlobal(wg, calka, mapa, blad);

    fclose(calka);
    fclose(mapa);
    fclose(blad);

    wg = 0.6;

    calka = fopen("calka0.6_glob.dat", "w");
    mapa = fopen("mapa0.6_glob.dat", "w");
    blad = fopen("blad0.6_glob.dat", "w");

    printf("Relaksacja globalna 0.6 iteracje");
    relGlobal(wg, calka, mapa, blad);

    fclose(calka);
    fclose(mapa);
    fclose(blad);

// zad 2

    double wl;

    wl = 1.0;
    calka = fopen("calka1.0_loc.dat", "w");
    printf("Relaksacja loaklna 1.0 iteracje");
    relLocal(wl, calka);
    fclose(calka);

    wl = 1.4;
    calka = fopen("calka1.4_loc.dat", "w");
    printf("Relaksacja loaklna 1.4 iteracje");
    relLocal(wl, calka);
    fclose(calka);

    wl = 1.8;
    calka = fopen("calka1.8_loc.dat", "w");
    printf("Relaksacja loaklna 1.8 iteracje");
    relLocal(wl, calka);
    fclose(calka);

    wl = 1.9;
    calka = fopen("calka1.9_loc.dat", "w");
    printf("Relaksacja loaklna 1.9 iteracje");
    relLocal(wl, calka);
    fclose(calka);

    return 0;
}