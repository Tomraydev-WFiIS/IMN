#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "mgmres.h" //załączenie nagłówkowego wyrzuca błąd o niezdefiniowanej referencji

#define delta 0.1

int N(int nx, int ny){
    return (nx+1)*(ny+1);
}

int j(int nx, int l){
    return (l/(nx+1));
}

int i(int nx, int l){
    return l-j(nx,l)*(nx+1);
}

double el(int nx, int l, int eps1, int eps2){
    if(i(nx,l) <= nx/2)
        return eps1;
    else
        return eps2;
}

double ro1(double x, double y, double xmax, double ymax, double sig) {
	return exp(-1*pow(x - 0.25*xmax, 2)/pow(sig, 2) - pow(y - 0.5*ymax, 2)/pow(sig, 2));
}


double ro2(double x, double y, double xmax, double ymax, double sig) {
	return -1*exp(-1*pow(x - 0.75*xmax, 2)/pow(sig, 2) - pow(y - 0.5*ymax, 2)/pow(sig, 2));
}

int dirichlet(double V1, double V2, double V3, double V4, int nx, int ny, int eps1, int eps2,
    int *ja, int *ia, double *a, double *b,double xmax, double ymax, FILE * macierz, FILE * wektor) {
//numeruje niezerowe
	int k = -1;
//liczba niezerowych el
	int nz_num = 0;

	for(int l = 0; l < N(nx, ny); ++l) {
		int brzeg = 0;  //0 srodek, 1 brzeg
		double vb = 0;  //potencjal na brzegu
		if(i(nx, l) == 0){
			brzeg = 1;
			vb = V1;
		}

		else if(i(nx, l) == nx){
			brzeg = 1;
			vb = V3;
		}

		else if(j(nx, l) == ny){
			brzeg = 1;
			vb = V2;
		}

		else if(j(nx, l) == 0){
			brzeg = 1;
			vb = V4;
		}

        b[l] = (-1)*(ro1(delta*i(nx,l), delta*j(nx,l),xmax,ymax,xmax/10) + ro2(delta*i(nx,l), delta*j(nx,l),xmax,ymax,xmax/10)); //sigma


		if(brzeg == 1)
			b[l] = vb;


		ia[l] = -1;

		if(l - nx - 1 > 0 && brzeg == 0){
			k++;
			if(ia[l] < 0)
                ia[l] = k;

			a[k] = el(nx,l,eps1,eps2)/(delta*delta);
			ja[k] = l - nx - 1;
		}
//poddiag
		if(l-1 > 0 && brzeg == 0) {
			k++;
			if(ia[l] < 0)
                ia[l] = k;
			a[k] = el(nx,l,eps1,eps2)/(delta*delta);
			ja[k] = l-1;
		}
//diag
		k++;
		if(ia[l] < 0)
            ia[l] = k;

		if(brzeg == 0)
			a[k] = -(2*el(nx,l,eps1,eps2)+el(nx,l+1,eps1,eps2)+el(nx,l+nx+1,eps1,eps2))/(delta*delta);
		else
			a[k] = 1;

		ja[k] = l;
//naddiag
		if(l < N(nx, ny) && brzeg == 0){
			k++;
			a[k] = el(nx,l+1,eps1,eps2)/(delta*delta);
			ja[k] = l + 1;
		}
//prawa skrajna przek
		if(l < N(nx, ny)-nx-1 && brzeg == 0){
			k++;
			a[k] = el(nx,l+nx+1,eps1,eps2)/(delta*delta);
			ja[k] = l + nx + 1;
		}
        if(l%5 == 0 && l != 0)
            fprintf(wektor, "\n");
        fprintf(wektor,"%d %d %d %f \n", l, i(nx,l), j(nx,l), b[l]);

    }

	nz_num = k+1;
	ia[N(nx, ny)] = nz_num;
    for(int z = 0; z < 5*N(nx,ny); z++)
        fprintf(macierz,"%d %0.f \n", z, a[z]);

    return nz_num;
}

void poisson(int nx, int ny, int eps1, int eps2, int V1, int V2, int V3, int V4, double xmax, double ymax, FILE * macierz, FILE * wektor, FILE * mapa = NULL){

//wektory
    double a[5*N(nx,ny)];
    int ja[5*N(nx,ny)];
    int ia[N(nx,ny)+1];
    double b[N(nx,ny)];
    double V[N(nx,ny)];
	//inicjalizacja a

    int nz_num = dirichlet(V1,V2,V3,V4,nx,ny,eps1,eps2,ja,ia,a,b,xmax,ymax, macierz, wektor);

    int itr_max = 500;
    int mr = 500;
    double tol_abs = pow(10,-8);
    double tol_rel = pow(10,-8);
    pmgmres_ilu_cr(N(nx,ny),nz_num,ia,ja,a,V,b,itr_max, mr,tol_abs,tol_rel);

//zapis mapy
	if(mapa){
		double tmp =0.;
		for(int z = 0; z < N(nx, ny); ++z){
			if(delta*j(nx,z) > tmp)
				fprintf(mapa,"\n");
			fprintf(mapa,"%f %f %f \n", delta*j(nx,z), delta*i(nx,z), V[z]);
			tmp = delta*j(nx,z);
		}
	}
}

int main(){
//definiowanie
    int nx = 4;
    int ny = 4;
    int eps1 = 1;
    int eps2 = 1;
    int V1 = 10;
    int V2 = -10;
    int V3 = 10;
    int V4 = -10;
    double xmax = 0;
    double ymax = 0;

//test
    FILE * macierz = fopen("macierz_test.dat", "w");
    FILE * wektor = fopen("wektor_test.dat", "w");
    poisson(nx,ny,eps1,eps2,V1,V2,V3,V4,xmax,ymax,macierz,wektor);
    fclose(macierz);
    fclose(wektor);
//zad 1
	nx = 50;
	ny = 50;
 	macierz = fopen("tmp.dat", "w");	//plik tmp dla wygody
    wektor = fopen("tmp.dat", "w");
	FILE * mapa = fopen("mapa50.dat", "w");
    poisson(nx,ny,eps1,eps2,V1,V2,V3,V4,xmax,ymax,macierz,wektor, mapa);
    fclose(mapa);

	nx = 100;
	ny = 100;

	mapa = fopen("mapa100.dat", "w");
    poisson(nx,ny,eps1,eps2,V1,V2,V3,V4,xmax,ymax,macierz,wektor, mapa);
    fclose(mapa);

	nx = 200;
	ny = 200;

	mapa = fopen("mapa200.dat", "w");
    poisson(nx,ny,eps1,eps2,V1,V2,V3,V4,xmax,ymax,macierz,wektor, mapa);
    fclose(mapa);
//zad2
	nx = 100;
	ny = 100;
	V1=V2=V3=V4=0;
	xmax = delta*nx;
	ymax = delta*ny;
//a) eps1 = eps2 = 1

	mapa = fopen("mapa2a.dat", "w");
    poisson(nx,ny,eps1,eps2,V1,V2,V3,V4,xmax,ymax,macierz,wektor, mapa);
    fclose(mapa);

//b) eps1 = 1 eps2 = 2
	eps1 = 1;
	eps2 = 2;

	mapa = fopen("mapa2b.dat", "w");
    poisson(nx,ny,eps1,eps2,V1,V2,V3,V4,xmax,ymax,macierz,wektor, mapa);
    fclose(mapa);

//c) eps1 = 1 eps2 = 10
	eps1 = 1;
	eps2 = 10;
	mapa = fopen("mapa2c.dat", "w");
    poisson(nx,ny,eps1,eps2,V1,V2,V3,V4,xmax,ymax,macierz,wektor, mapa);
    fclose(mapa);

    fclose(macierz);
    fclose(wektor);

}