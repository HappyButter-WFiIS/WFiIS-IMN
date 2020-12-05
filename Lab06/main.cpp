#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "mgmres.h"


double getRho1(double x, double xMax, double y, double yMax, double sigma){
    return 1.0 * exp(-pow(x-0.25*xMax, 2)/pow(sigma,2) - pow(y-0.5*yMax, 2)/pow(sigma,2));
}

double getRho2(double x, double xMax, double y, double yMax, double sigma){
    return -1.0 * exp(-pow(x-0.75*xMax, 2)/pow(sigma,2) - pow(y-0.5*yMax, 2)/pow(sigma,2));
}

double getRhoSum(double x, double xMax, double y, double yMax, double sigma){
    return getRho1(x, xMax, y, yMax, sigma) + getRho2(x, xMax, y, yMax, sigma);
}

int getJ(int l, int nx){
    return floor(l/(nx+1));
}

int getI(int l, int nx){
    return l-getJ(l,nx)*(nx+1);
}

double getEpsilonl(int i, int nx, int epsilon1, int epsilon2){
    if(i<=nx/2)
        return epsilon1;
    return epsilon2;
}

int getN(int nx, int ny){
    return (nx+1)*(ny+1);
}


void saveVectorBData(int N, int nx, double b[], FILE * vectorB){
    for(int m = 0; m < N; m++){
        double i = getI(m, nx);
        double j = getJ(m, nx);
        fprintf(vectorB,"%d %.1f %.1f %f \n", m, i, j, b[m]);
    }
}

void saveMatrixAData(int N, int nx, double a[], FILE * matrixA){
    int matrixSize = 5*N;

    for(int m = 0; m < matrixSize; m++){
        double i = getI(m, nx);
        double j = getJ(m, nx);
        fprintf(matrixA,"%d %.1f %.1f %0.f \n", m, i, j, a[m]);
    }
}

int Dirichlet(  int nx, int ny, int epsilon1, int epsilon2,
                int V1, int V2, int V3, int V4,
                double xMax, double yMax, double sigma, double delta,
                double a[], double b[], int ia[], int ja[], 
                int N, FILE* matrixA = nullptr, FILE* vectorB = nullptr){
    
    int k = -1;

    for(int it = 0; it<N; it++){
        int brzeg = 0;
        double vb = 0.0;

        if(getI(it,nx) == 0){
            brzeg = 1;
            vb = V1;
        }

        if(getJ(it,nx) == ny){
            brzeg = 1;
            vb = V2;
        }

        if(getI(it,nx) == nx){
            brzeg = 1;
            vb = V3;
        }

        if(getJ(it,nx) == 0){
            brzeg = 1;
            vb = V4;
        }

        b[it] = -(getRhoSum(delta*getI(it,nx), xMax, delta*getJ(it,nx),yMax, sigma));

        if(brzeg == 1){
            b[it] = vb;
        }

        ia[it] = -1;

        // lewa skrajna przekątna
        if(it-nx-1 >= 0 && brzeg == 0){
            k++;
            if(ia[it]<0)
                ia[it] = k;
            a[k] = getEpsilonl(getI(it,nx), nx, epsilon1, epsilon2) / pow(delta,2);
            ja[k] = it - nx - 1;
        }

        // poddiagonala
        if(it-1 >= 0 && brzeg == 0){
            k++;
            if(ia[it]<0)
                ia[it] = k;
            a[k] = getEpsilonl(getI(it,nx), nx, epsilon1, epsilon2) / pow(delta,2);
            ja[k] = it - 1;
        } 

        // diagonala
        k++;
        if(ia[it]<0)
            ia[it] = k;

        if(brzeg == 0){
            a[k] = -1.0 * (2.0 * getEpsilonl(getI(it,nx), nx, epsilon1, epsilon2) + 
                    getEpsilonl(getI(it+1,nx), nx, epsilon1, epsilon2) + 
                    getEpsilonl(getI(it+nx+1,nx), nx, epsilon1, epsilon2) ) 
                    / pow(delta,2);    
        }else{
            a[k] = 1;
        }
        ja[k] = it;
        
        // naddiagonala
        if(it < N && brzeg == 0){
            k++;
            a[k] = getEpsilonl(getI(it+1,nx), nx, epsilon1, epsilon2) / pow(delta,2);
            ja[k] = it + 1;
        }

        // prawa skrajna przekątna
        if(it < N-nx-1 && brzeg == 0){
            k++;
            a[k] = getEpsilonl(getI(it+nx+1,nx), nx, epsilon1, epsilon2) / pow(delta,2);
            ja[k] = it + nx + 1;
        }
    }


    int nz_num = k+1;
	ia[N] = nz_num;


    // zapisz vector b jeśli plik nie jest nullptr
    if(vectorB){
        saveVectorBData(N, nx, b, vectorB);
    }

    // zapisz macierz A jeśli plik nie jest nullptr
    if(matrixA){
        saveMatrixAData(N,nx,a,matrixA);
    }

    return nz_num;
}

void saveMapData(int N, int nx, double delta, double V[], FILE * mapData){
    double makeBreak = 0.0;
    for(int m = 0; m < N; m++){
        double i = getI(m, nx);
        double j = getJ(m, nx);

        if(delta*j > makeBreak)
            fprintf(mapData,"\n");

        fprintf(mapData,"%f %f %f \n", delta*i, delta*j, V[m]);

        makeBreak = delta*j;
    }
}


void calculate( int nx, int ny, int epsilon1, int epsilon2,
                int V1, int V2, int V3, int V4,
                double xMax, double yMax, double sigma, double delta,
                FILE * matrixA = nullptr, FILE * vectorB = nullptr, FILE * mapData = nullptr){


    // local init conditions
    int itr_max = 500;
    int mr = 500;
    double tol_abs = 1e-8;
    double tol_rel = 1e-8;

    int N = getN(nx,ny);
    
    // create vectors    
    double a[5*N];
    int ja[5*N];
    int ia[N+1];
    double b[N];
    double V[N];

    // fill vectors
    int nz_num = Dirichlet(nx, ny, epsilon1, epsilon2,
                V1, V2, V3, V4,
                xMax, yMax, sigma, delta,
                a, b, ia, ja, 
                N, matrixA, vectorB);
    // solve
    pmgmres_ilu_cr(N,nz_num,ia,ja,a,V,b,itr_max, mr,tol_abs,tol_rel);

    // save result if needed
	if(mapData){
        saveMapData(N, nx, delta, V, mapData);
	}
}


int main(){

    // init conditions
    double delta = 0.1;
    
    int nx = 4;
    int ny = 4;
    int epsilon1 = 1;
    int epsilon2 = 1;
    int V1 = 10;
    int V3 = 10;
    int V2 = -10;
    int V4 = -10;

    double xMax = 0.0;
    double yMax = 0.0;
    double sigma = xMax/10.;

    // Zadanie nr 1

    FILE * matrixA = fopen("matrixA.dat", "w+");
    FILE * vectorB = fopen("vectorB.dat", "w+");

    fprintf(matrixA,"#l i j a[l] \n\n");
    fprintf(vectorB,"#l i j b[l] \n\n");

    calculate(  nx, ny, epsilon1, epsilon2,
                V1, V2, V3, V4,
                xMax, yMax, sigma, delta,
                matrixA, vectorB);

    fclose(matrixA);
    fclose(vectorB);
    
    // Zadanie nr 2
    // a)
    nx = 50;
	ny = 50;
    
	FILE * mapData = fopen("zad2_50.dat", "w+");
    calculate(  nx, ny, epsilon1, epsilon2,
                V1, V2, V3, V4,
                xMax, yMax, sigma, delta,
                nullptr, nullptr, mapData);
    fclose(mapData);

    // b)
    nx = 100;
	ny = 100;

	mapData = fopen("zad2_100.dat", "w+");
    calculate(  nx, ny, epsilon1, epsilon2,
                V1, V2, V3, V4,
                xMax, yMax, sigma, delta,
                nullptr, nullptr, mapData);
    fclose(mapData);

    // c)
	nx = 200;
	ny = 200;

	mapData = fopen("zad2_200.dat", "w+");
    calculate(  nx, ny, epsilon1, epsilon2,
                V1, V2, V3, V4,
                xMax, yMax, sigma, delta,
                nullptr, nullptr, mapData);
    fclose(mapData);

    // Zadanie nr 3
    // new init conditions
    V1 = 0;
    V3 = 0;
    V2 = 0;
    V4 = 0;

    nx = 100;
	ny = 100;
    xMax = delta*((double)nx);
	yMax = delta*((double)ny);
    sigma = xMax/10.;

    // a)
    epsilon1 = 1;
    epsilon2 = 1;

	mapData = fopen("zad3_a.dat", "w+");
    calculate(  nx, ny, epsilon1, epsilon2,
                V1, V2, V3, V4,
                xMax, yMax, sigma, delta,
                nullptr, nullptr, mapData);
    fclose(mapData);

    // b)
    epsilon1 = 1;
    epsilon2 = 2;

    mapData = fopen("zad3_b.dat", "w+");
    calculate(  nx, ny, epsilon1, epsilon2,
                V1, V2, V3, V4,
                xMax, yMax, sigma, delta,
                nullptr, nullptr, mapData);
    fclose(mapData);

    // c)
    epsilon1 = 1;
    epsilon2 = 10;
    
    mapData = fopen("zad3_c.dat", "w+");
    calculate(  nx, ny, epsilon1, epsilon2,
                V1, V2, V3, V4,
                xMax, yMax, sigma, delta,
                nullptr, nullptr, mapData);
    fclose(mapData);
}