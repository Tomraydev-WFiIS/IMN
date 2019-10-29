#include <cstdio>
#include <cstdlib>
#include <cmath>

//predefiniowane stale podane w poleceniu - we wszystkich zadaniach takie same
#define beta 0.001
#define N 500
#define gamma 0.1
#define t_max 100
#define dt 0.1
#define mi_max 20

//fun niezalezna od t
double f(double u){
    double alfa = beta*N - gamma;
    return (alfa*u - beta*u*u);
}

//pochodna powyzszej funkcji
double df(double u){
    double alfa = beta*N - gamma;
    return (alfa- 2*beta*u);
}

//zadanie 1
void metoda_picarda(FILE * file, double u0, double TOL){

    double t = 0.;
    double un =u0;
    double un1 = u0;
    double umi = 0.;
    int mi_it = 0;
    //maksymalny nr iteracji w petli wewnetrznej
    int max_iter = 0;

    //w chwili 0.0 zakladam, ze jest jeden zarazony
    fprintf(file, "%.1f %f %f \n" , t, un1, N-un1);


    for(t = 0.1 ; t < t_max; t+=dt ){

        mi_it = 0;
        un = un1;
        umi = 0.;   //warunek w while nie bedzie spelniony wiec zawsze wejdzie do petli.
        // Warunek startowy ustanawiam w petli wewnetrznej (umi = un)

        while(fabs(un1 - umi) > TOL && mi_it++ < mi_max){
            umi = un1;  //w pierwszej iteracji un1 = un - warunek startowy
            un1 = un + (dt/2.)*(f(un)+f(umi));
             if(mi_it > max_iter)
                max_iter = mi_it;
        }

        fprintf(file, "%.1f %f %f \n" , t, un1, N-un1);
    }

    printf("metoda Picarda  %d \n", max_iter);

}

//zadanie 2
void iteracja_newtona(FILE * file, double u0, double TOL){

    double t = 0.;
    double un =u0;
    double un1 = u0;
    double umi = 0.;
    int mi_it = 0;
    //maksymalna nr iteracji w petli wewnetrznej
    int max_iter = 0;

    //w chwili 0.0 zakladam, ze jest jeden zarazony
    fprintf(file, "%.1f %f %f \n" , t, un1, N-un1);

    for(t = 0.1  ; t < t_max; t+=dt ){

        mi_it = 0;
        un = un1;
        umi = 0.; //analogicznie do poprzedniego zadania

        while(fabs(un1 - umi) > TOL && mi_it++ < mi_max  ){
            umi = un1; //w pierwszej iteracji un1 = un - warunek startowy
            un1 = umi - (umi - un - (dt/2.)*(f(un)+f(umi)))/(1 - (dt/2.)*df(umi));
            if(mi_it > max_iter)
                max_iter = mi_it;
        }
        fprintf(file, "%.1f %f %f \n" , t, un1, N-un1);
    }
    printf("iteracja newtona  %d \n", max_iter);
}


//zadanie 3
double F1(double U1, double U2, double un,double a11, double a12){
    double alfa = beta*N - gamma;
    return (U1 - un - dt*(a11*(alfa*U1-beta*U1*U1) + a12*(alfa*U2 - beta*U2*U2)));
}
double F2(double U1, double U2, double un,double a21, double a22){
    double alfa = beta*N - gamma;
    return (U2 - un - dt*(a21*(alfa*U1-beta*U1*U1) + a22*(alfa*U2 - beta*U2*U2)));
}
//elementy macierzowe, kolejne kombinacje pochodnych
double m11(double U1, double a11){
    double alfa = beta*N - gamma;
    return (1-dt*a11*(alfa - 2*beta*U1));
}

double m12(double U2, double a12){
    double alfa = beta*N - gamma;
    return ((-1)*dt*a12*(alfa - 2*beta*U2));
}

double m21(double U1, double a21){
    double alfa = beta*N - gamma;
    return ((-1)*dt*a21*(alfa - 2*beta*U1));
}

double m22(double U2, double a22){
    double alfa = beta*N - gamma;
    return (1-dt*a22*(alfa - 2*beta*U2));
}


void niejawna_metoda_RK2(FILE *file, double u0, double TOL){
    //f nie jest zalezna od t wiec wartosci c1 i c2 nie sa potrzebne


    double a11 = 0.25;
    double a12 = 0.25 - sqrt(3)/6;

    double a21 = 0.25 + sqrt(3)/6;
    double a22 = a11;

    double b = 0.5; //b1 = b2

    double U1 = 0., U2 = 0., U1mi = 0., U2mi = 0.;

    double t = 0.;
    double un =u0;
    double un1 = u0;
    int mi_it = 0;
    //maksymalna nr iteracji w petli wewnetrznej
    int max_iter = 0;

    //w chwili 0.0 zakladam, ze jest jeden zarazony
    fprintf(file, "%.1f %f %f \n" , t, un1, N-un1);

    for(t = 0.1  ; t < t_max; t+=dt ){

        mi_it = 0;
        un = un1;
        U1mi = 0.;
        U2mi = 0.;
        U1 = un;
        U2 = un;

        while((fabs(U1 - U1mi) > TOL || fabs(U2 - U2mi) > TOL) && mi_it++ < mi_max  ){
            U1mi = U1;
            U2mi = U2;

            double dtU1 = (F2(U1,U2,un,a21,a22)*m12(U2,a12)-F1(U1,U2,un,a11,a12)*m22(U2,a22))/
                                                (m11(U1,a11)*m22(U2,a22)-m12(U2,a12)*m21(U1,a21));

            double dtU2 = (F1(U1,U2,un,a11,a12)*m21(U1,a21)-F2(U1,U2,un,a21,a22)*m11(U1,a11))/
                                                (m11(U1,a11)*m22(U2,a22)-m12(U2,a12*m21(U1,a21)));

            U1 = U1mi+dtU1;
            U2 = U2mi + dtU2;

            if(mi_it > max_iter)
                max_iter = mi_it;
        }
        un1 = un + dt*(b*f(U1) + b*f(U2));
        fprintf(file, "%.1f %f %f \n" , t, un1, N-un1);

    }
    printf("iteracja RK2  %d \n", max_iter);

}

int main(){

    double  u0 = 1.;
    double TOL = pow(10,-6);
    FILE * file;

//zadanie 1
    file = fopen("zad1_rozw.dat", "w");
    metoda_picarda(file, u0, TOL);
    fclose(file);

//zadanie 2
    file = fopen("zad2_rozw.dat", "w");
    iteracja_newtona(file, u0, TOL);
    fclose(file);

//zad3
    file = fopen("zad3_rozw.dat", "w");
    niejawna_metoda_RK2(file,u0,TOL);
    fclose(file);

}