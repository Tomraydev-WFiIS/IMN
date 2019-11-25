#include <stdlib.h>
#include <cmath>
#include <ctime>    
#include <iostream>


#define delta 0.2
#define nx 128
#define ny 128
#define x_max (delta)*(nx)
#define y_max (delta)*(ny)
#define TOL 1e-8    // tolerance

#define VB1(y) (sin(M_PI * (y)/(y_max) ))   // left boundary condition
#define VB2(x) (-sin(2*M_PI * (x)/(x_max) ))   // top boundary condition
#define VB3(y) (sin(M_PI * (y)/(y_max) ))   // right boundary condition
#define VB4(x) (sin(2*M_PI * (x)/(x_max) ))   // bottom boundary condition


double rho[nx+1][ny+1]; // charge density map

double rho1 (int i, int j){
    return exp(
        -pow((i*delta - 0.35*x_max),2)/pow(sigma_x,2)
        -pow((j*delta - 0.5*y_max),2)/pow(sigma_y,2)
    );
}

double rho2 (int i, int j){
    return -exp(
        -pow((i*delta - 0.65*x_max),2)/pow(sigma_x,2)
        -pow((j*delta - 0.5*y_max),2)/pow(sigma_y,2)
    );
}


int main(void){

}