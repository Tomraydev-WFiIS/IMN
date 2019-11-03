#include <cstdio>
#include <cstdlib>
#include <cmath>

//parametry startowe predefiniowane
#define x0 0.01
#define v0 0.0
#define dt0 1.0
#define S 0.75
#define p 2.0

//pomocnicze do algorytmu

//nie trzeba się bawic z pamiecia i zwracaniem tablicy ( można też std::pair)
struct pair{
    double x;
    double v;
};

// można też std::max
double max(double x, double y){
    return x > y ? x : y;
}

//do trapezow i RK2

//niezalezne od xn
double f(double vn){
    return vn;
}

double g(double xn, double vn, int alfa){
    return (alfa*(1-xn*xn)*vn - xn);
}

//metoda RK2 - zadanie 1
pair m_RK2(double xn, double vn, double dt, int alfa){
    double k1x = f(vn);
    double k1v = g(xn, vn, alfa);

    double k2x = f(vn+dt*k1v);
    double k2v = g(xn + dt*k1x, vn + dt*k1v, alfa);

    pair toRet;
    toRet.x = xn + (dt/2)*(k1x + k2x);
    toRet.v = vn + (dt/2)*(k1v + k2v);
    return toRet;

}

//metoda trapezow - zadanie 2

double F(double xn1, double xn, double vn1, double vn, double dt){

    return xn1 - xn - (dt/2)*(f(vn)+f(vn1));
}

double G(double xn1, double xn, double vn1, double vn, double dt,int alfa){

    return vn1 - vn - (dt/2)*(g(xn,vn, alfa)+g(xn1,vn1, alfa));
}

//elementy macierzowe
double a11(){
    return 1.;
}

double a12(double dt){
    return (-1)*(dt/2);
}

double a21(double dt, double xn1k, double vn1k, double alfa){
    return (-1)*(dt/2)*((-2)*alfa*xn1k*vn1k-1);
}

double a22(double dt, double xn1k, double alfa){
    return (1 - (dt/2)*alfa*(1-xn1k*xn1k));
}


//zadanie 2
pair m_trapezow(double xn, double vn, double dt, int alfa){
    double delta = pow(10,-10);

    double xn1k = xn;
    double vn1k = vn;
    double xn1 = xn;
    double vn1 = vn;
    double dx, dv;

    do{
        dx = 0.;
        dv = 0.;
        vn1k = vn1;
        xn1k = xn1;

        dx = ((-1)*F(xn1,xn, vn1, vn, dt)*a22(dt,xn1k,alfa)-(-1)*G(xn1,xn,vn1,vn,dt, alfa)*a12(dt))/
                                            (a11()*a22(dt,xn1k,alfa) - a12(dt)*a21(dt,xn1k,vn1k,alfa));

        dv = (a11()*(-1)*G(xn1,xn,vn1,vn,dt, alfa) - a21(dt,xn1k,vn1k,alfa)*(-1)*F(xn1,xn,vn1,vn,dt))/
                                            (a11()*a22(dt,xn1k,alfa) - a12(dt)*a21(dt,xn1k,vn1k,alfa));

        xn1 = xn1k + dx;
        vn1 = vn1k + dv;
    }while(fabs(dx) > delta || fabs(dv) > delta);


    pair toRet;

    toRet.x = xn + (dt/2)*(f(vn) + f(vn1));
    toRet.v = vn + (dt/2)*(g(xn,vn, alfa) + g(xn1,vn1, alfa));

    return toRet;
}

//algorytm uniwersalny, przyjmuje wskaznik na schemat numeryczny
void algorytm_kontroli_kroku(FILE * file, double TOL, pair (*schemat_numeryczny)(double xn, double vn, double dt, int alfa)){
//mamy przekazac alfe jako argument, dlatego nie predefiniowana
    int alfa = 5;

    double t = 0.0;
    double  dt = dt0;
    double xn = x0;
    double vn = v0;
    double tmax = 40.0;

    double Ex = 0.;
    double Ev = 0.;
    pair tmp2;
    pair tmp1;

    //dla t = 0.0
    fprintf(file, "%f %f %f %f \n", t, dt, xn, vn);
    do{
        //stawiamy dwa kroki dt
        tmp2 = schemat_numeryczny(xn,vn,dt,alfa);
        tmp2 = schemat_numeryczny(tmp2.x,tmp2.v,dt,alfa);

        //stawiamy jeden krok 2dt
        tmp1 = schemat_numeryczny(xn,vn,2*dt,alfa);

        //liczymy Ex Ev
            Ex = (tmp2.x - tmp1.x)/(pow(2,p) - 1);
            Ev = (tmp2.v - tmp1.v)/(pow(2,p) - 1);

        if(max(fabs(Ex),fabs(Ev)) < TOL ){
            t =t+2*dt;
            xn = tmp2.x;
            vn = tmp2.v;
            fprintf(file, "%f %f %f %f \n", t, dt, xn, vn);

        }
        dt = (pow(S*TOL/(max(fabs(Ex),fabs(Ev))),(1/(p+1)))*dt);

    }while(t < (tmax - dt));

}


int main(){

    FILE * file;

    double TOL = pow(10,-2);

// zadanie 1
    file = fopen("RK2.dat", "w");
    algorytm_kontroli_kroku(file, TOL, m_RK2);

    fprintf(file, "\n\n");

    TOL = pow(10,-5);
    algorytm_kontroli_kroku(file, TOL, m_RK2);

    fclose(file);


// zadanie 2

    file = fopen("trapezy.dat", "w");
    TOL = pow(10,-2);
    algorytm_kontroli_kroku(file, TOL, m_trapezow);

    fprintf(file, "\n\n");

    TOL = pow(10,-5);
    algorytm_kontroli_kroku(file, TOL, m_trapezow);

    fclose(file);
}