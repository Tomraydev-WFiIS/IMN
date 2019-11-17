#include <stdlib.h>
#include <cmath>
#include <ctime>    
#include <iostream>

// 1. Przyjmujemy wartości parametrów: ∆ = 0.2, nx = 128, ny = 128, xmax = ∆ · nx, ymax = ∆ · ny,
// T OL = 10−8 oraz warunki brzegowe Dirichleta:

#define delta 0.2
#define nx 128
#define ny 128
#define x_max (delta)*(nx)
#define y_max (delta)*(ny)
#define TOL 1e-8    // tolerance

#define VB1(y) (sin(M_PI * (y)/(y_max) ))   // left boundary condition
#define VB2(y) (-sin(2*M_PI * (x)/(x_max) ))   // top boundary condition
#define VB3(y) (sin(M_PI * (y)/(y_max) ))   // right boundary condition
#define VB4(y) (sin(2*M_PI * (x)/(x_max) ))   // bottom boundary condition

int main(void){
    // 2. Rozwiązać równanie Poissona z zadanymi WB metodą wielosiatkową dla k = 16, 8, 4, 2, 1. Dla
    // każdego k po spełnieniu warunku stopu sporządzić mapę potencjału (5 map). (60 pkt) Dla
    // każdego k zapisać do pliku wartości całki funkcjonalnej w funkcji numeru iteracji. Sporządzić
    // wykres zmian S^(k)(it) dla wszystkich k na jednym rysunku. (40 pkt)
    FILE * f_integral = fopen("integral.dat", "w");
    FILE * f_maps = fopen("maps.dat", "w");


    fclose(f_integral);
    fclose(f_maps);

    // Uwaga 1: Wszystkie obliczenia wykonujemy korzystając z jednej tablicy potencjału (jak dla najgęstszej
    // siatki), w której poruszamy się z aktualnym krokiem k.

    // Uwaga 2: Warunki brzegowe wyznaczamy tylko raz - przed rozpoczęciem relaksacji na najrzadszej
    // siatce. WB określamy dla każdego węzła brzegowego (jak dla k = 1). Po określeniu WB, zerujemy
    // potencjał w każdym węźle (k = 1) wewnątrz obszaru (start metody).

    // Uwaga 3: Po uzyskaniu samouzgodnienia na siatce o indeksie k, zagęszczamy siatkę tj. w nowych
    // węzłach (czerwonych) wpisujemy wartości interpolowane. Jest to potencjał startowy dla relaksacji na
    // gęstszej siatce, gdyż stanowi on (na ogół) dobre przybliżenie dokładnego rozwiązania.
        

    return 0;
}






