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
#define VB2(x) (-sin(2*M_PI * (x)/(x_max) ))   // top boundary condition
#define VB3(y) (sin(M_PI * (y)/(y_max) ))   // right boundary condition
#define VB4(x) (sin(2*M_PI * (x)/(x_max) ))   // bottom boundary condition

double sum(double V[][ny + 1], int k){
    double S = 0.0;
	for(int i = 0; i <= nx-k; i += k){
		for(int j = 0; j <= ny-k; j += k){
			S += pow(k*delta, 2)/2.0*
                (
                    pow( (V[i+k][j] - V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i][j+k])/(2*k*delta), 2) +
					pow( (V[i][j+k] - V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i+k][j])/(2*k*delta), 2)
                );
		}
	}
	return S;
}


void multigrid_relaxation(FILE * f_maps, FILE * f_integral){
    std::clock_t start;
    double duration = 0;
    double total_duration;
    int it = 0;
    int total_it = 0;
    double V[nx+1][ny+1] = {0.0}; // map of potential

    // boundary conditions
    for(int i = 0; i <= nx; i++){
        V[i][ny] = VB2(i*delta); // top
        V[i][0] = VB4(i*delta); // bottom
    }
    for(int i = 0; i <= ny; i++){
        V[0][i] = VB1(i*delta); // left
        V[nx][i] = VB3(i*delta); // right
    }

    int k = 16;
    double S1 = sum(V, k);
    double S = S1;
    double diff;
    while(k > 0){
        start = std::clock();
        it = 0;
        // relax the net
        do{
            for(int i = k; i <= nx-k; i+=k){
                for(int j = k; j <= ny-k; j+=k){
                    V[i][j] = 0.25*(V[i+k][j] + V[i-k][j] + + V[i][j+k] + + V[i][j-k]);
                }
            }

            // Calculate the stop condition.
            S = S1;
            S1 = sum(V, k);
            fprintf(f_integral, "%d %d %lf \n",k, total_it + it, S1);
            diff = (fabs(S1 - S) / S);
        } while ( diff > TOL && ++it );
        fprintf(f_integral, "\n\n");
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        printf("k = %d\titerations = %d\tduration=%lf s\n", k, it, duration);

        total_it += it;
        total_duration += duration;


        // save the map
        for(int i = 0; i <= nx; i += k){
            for(int j = 0; j <= ny; j += k){
                fprintf(f_maps, "%lf %lf %lf \n", delta*i, delta*j, V[i][j]);
            }
            fprintf(f_maps, "\n");
        }
        fprintf(f_maps, "\n");


        // finer net
        for (int i = 0; i <= nx-k; i += k){
            for (int j = 0; j <= ny-k; j += k){
                                V[i+k/2][j+k/2] = 0.25*(V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]); // center
                if(i != nx-k)   V[i+k][j+k/2] = 0.5*(V[i+k][j] + V[i+k][j+k]); // right
                if(j != ny-k)   V[i+k/2][j+k] = 0.5*(V[i][j+k] + V[i+k][j+k]); // top
                if(j != 0)      V[i+k/2][j] = 0.5*(V[i][j] + V[i+k][j]); // bottom
                if(i != 0)      V[i][j+ k/2] = 0.5*(V[i][j] + V[i][j+k]); // left
            }
        }
        k /= 2;
    } // endwhile
    printf("total_iterations = %d\ttotal_duration=%lf s\n", total_it, total_duration);
}

int main(void){
    // 2. Rozwiązać równanie Poissona z zadanymi WB metodą wielosiatkową dla k = 16, 8, 4, 2, 1. Dla
    // każdego k po spełnieniu warunku stopu sporządzić mapę potencjału (5 map). (60 pkt) Dla
    // każdego k zapisać do pliku wartości całki funkcjonalnej w funkcji numeru iteracji. Sporządzić
    // wykres zmian S^(k)(it) dla wszystkich k na jednym rysunku. (40 pkt)
    FILE * f_maps = fopen("maps.dat", "w");
    FILE * f_integral = fopen("integral.dat", "w");
    multigrid_relaxation(f_maps, f_integral);
    fclose(f_maps);
    fclose(f_integral);
    return 0;
}
    // Uwaga 1: Wszystkie obliczenia wykonujemy korzystając z jednej tablicy potencjału (jak dla najgęstszej
    // siatki), w której poruszamy się z aktualnym krokiem k.

    // Uwaga 2: Warunki brzegowe wyznaczamy tylko raz - przed rozpoczęciem relaksacji na najrzadszej
    // siatce. WB określamy dla każdego węzła brzegowego (jak dla k = 1). Po określeniu WB, zerujemy
    // potencjał w każdym węźle (k = 1) wewnątrz obszaru (start metody).

    // Uwaga 3: Po uzyskaniu samouzgodnienia na siatce o indeksie k, zagęszczamy siatkę tj. w nowych
    // węzłach (czerwonych) wpisujemy wartości interpolowane. Jest to potencjał startowy dla relaksacji na
    // gęstszej siatce, gdyż stanowi on (na ogół) dobre przybliżenie dokładnego rozwiązania.
