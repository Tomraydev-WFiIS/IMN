#include <cmath>
#include <iostream>
#include <string>

// 3. Przyjąć parametry startowe:
//  x0 = 0.01,
//  v0 = 0,
//  ∆t0 = 1,
//  S = 0.75,
//  p = 2 (rząd dokładności obu metod),
//  tmax = 40, α = 5.

#define x0 0.01
#define v0 0.0
#define dt0 1.0
#define S 0.75
#define p 2.0
double TMAX = 40;
double ALFA = 5;


// 1. Zaprogramować metody: trapezów i RK2 jako dwie osobne procedury.
//  Do procedur przekazujemy: xn, vn, ∆t, α; a mają zwracać: xn+1, vn+1.
std::pair<double,double> RK2(double xn, double vn, double dt, double alfa){
    double k1x = vn;
    double k1v = alfa*(1-xn*xn)*vn - xn;

    double k2x = vn + dt*k1v;
    double k2v = alfa * (1.0 - pow((xn + dt*k1x), 2.0)) * (vn + dt*k1v) - (xn + dt*k1x);

    double xn1 = xn + dt/2.0 * (k1x + k2x);
    double vn1 = vn + dt/2.0 * (k1v + k2v);
    return std::pair<double,double>(xn1, vn1);
}

double g(double x, double v){
    return ALFA*(1.0 - x*x)*v - x;
}

double a11(){
    return 1.0;
}
double a12(double dt){
    return -dt/2.0;
}
double a21(double dt, double xn1, double vn1){
    return (-dt/2.0)*(-2.0*ALFA*xn1*vn1 - 1.0);
}
double a22(double dt, double xn1){
    return 1.0 - (dt/2.0)*ALFA*(1.0 - xn1*xn1);
}

double F(double xn, double xn1,double vn, double vn1, double dt){
    return xn1 - xn - (dt/2.0)*(vn+vn1);
}
double G(double xn, double xn1,double vn, double vn1, double dt){
    return vn1 - vn - (dt/2.0)*(g(xn,vn) + g(xn1, vn1));
}


std::pair<double,double> trapezy(double xn, double vn, double dt, double alfa){
    double delta = 1e-10;
    double dx = 0.0;
    double dv = 0.0;
    double xn1 = xn;
    double vn1 = vn;
    // int i = 0;
    do{
        dx = (-F(xn,xn1,vn,vn1,dt)*a22(dt,xn1) - (-G(xn,xn1,vn,vn1,dt))*a12(dt) )
            / (a11() * a22(dt,xn1) - a12(dt) * a21(dt,xn1,vn1) );

        dv = (a11()*(-G(xn,xn1,vn,vn1,dt)) - a21(dt,xn1,vn1)*(-F(xn,xn1,vn,vn1,dt)) )
            / (a11() * a22(dt,xn1) - a12(dt) * a21(dt,xn1,vn1) );

        xn1 = xn1+dx;
        vn1 = vn1+dv;
        // printf("dx=%g\tdv=%g\n", dx, dv);
        // i++;
    } while (fabs(dx) > delta || fabs(dv) > delta);
    
    // printf("%d\txn1=%g\tvn1=%g\n", i, xn1, vn1);
    return std::pair<double,double>(xn1, vn1);
}

// 2. Zaimplementować algorytm kontroli kroku czasowego.
void step_control(FILE * f, double TOL, std::pair<double,double> (*next_step)(double xn, double vn, double dt, double alfa)){
    // initialize
    double t = 0.0;
    double dt = dt0;
    double xn = x0;
    double vn = v0;
    double tmax = TMAX;
    double Ex = 0.0;
    double Ev = 0.0;
    std::pair<double,double>two_steps(xn, vn);
    std::pair<double,double>one_step(xn, vn);

    do{
        // make two steps of length dt
        two_steps = next_step(xn, vn, dt, ALFA);
        two_steps = next_step(two_steps.first, two_steps.second, dt, ALFA);

        // make one step of length 2*dt
        one_step = next_step(xn, vn, 2*dt, ALFA);

        // calculate the errors
        Ex = (two_steps.first - one_step.first) / (pow(2.0, p) - 1.0);
        Ev = (two_steps.second - one_step.second) / (pow(2.0, p) - 1.0);

        if(std::max(fabs(Ex), fabs(Ev)) < TOL){
            t += 2*dt;
            xn = two_steps.first;
            vn = two_steps.second;
            fprintf(f, "%g\t%g\t%g\t%g\t\n", t, dt, xn, vn);
        }
        dt = pow(
            S * TOL / (std::max(fabs(Ex), fabs(Ev))),
            (1.0/(p+1.0))
        ) * dt;

    } while (t < tmax);
    return;
}


int main(void)
{
// 4. Rozwiązać równanie (1) metodą RK2 z kontrolą kroku czasowego dla parametru TOL = 10−2; 10−5.
// Wykonać rysunki: x(t), v(t), ∆t(t) i v(x).
// Wykresy tej samej wielkości dla obu wartości TOL umieścić na jednym rysunku (będą 4 rysunki). (50 pkt)
    FILE * f1 = fopen("RK2.dat", "w");
    double TOL = 1e-2;
    step_control(f1, TOL, RK2);
    fprintf(f1,"\n\n");

    TOL = 1e-5;
    step_control(f1, TOL, RK2);
    fclose(f1);

// 5. Rozwiązać równanie (1) metodą trapezów z kontrolą kroku czasowego dla parametru TOL = 10−2 ; 10−5 .
// Wykonać rysunki: x(t), v(t), ∆t(t) i v(x).
// Wykresy tej samej wielkości dla obu wartości TOL umieścić na jednym rysunku (będą 4 rysunki). (50 pkt)
    FILE * f2 = fopen("trapezy.dat", "w");
    TOL = 1e-2;
    step_control(f2, TOL, trapezy);
    fprintf(f2,"\n\n");

    TOL = 1e-5;
    step_control(f2, TOL, trapezy);
    fclose(f2);
}
