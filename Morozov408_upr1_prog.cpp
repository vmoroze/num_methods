#include <iostream>
#include <cmath>

double K_finder(double a, int r);
double S_finder(double a, double *k, double *w, double *d, int l, int r, int m);
double summF(double a, double b, double k, double *w, double *d, int r, int m);

int main(void){

    double *massW, *massD;
    int m = 5;
    double a;
    double *massK, *massDelta, *massS, *massRel;
    int num;
    int N = -1;
    int r, l;

    a = sqrt(3);
    num = 10;

    try{
		massW = new double[m];
        massD = new double[m];
        massK = new double[num];
        massS = new double[num];
        massDelta = new double[num];
        massRel = new double[4];
	} catch(...){
		std::cout << "some trouble with memory" << std::endl;
		return -1;
	}

    massW[0] = 7./90; massW[1] = 32./90; massW[2] = 12./90;
    massW[3] = 32./90; massW[4] = 7./90;
    massD[0] = -1.; massD[1] = -1./2; massD[2] = 0.;
    massD[3] = 1./2; massD[4] = 1.;

    for(r = 0; r < num; r++){
        massK[r] = K_finder(a, r);
    }
    for(r = 0; r < num; r++){
        massS[r] = S_finder(a, massK, massW, massD, 0, r, m);
    }
    for(r = 0; r < num; r++){
        massDelta[r] = fabs(1 - massS[r]);
    }

    printf("\n");
    printf("r\tdelta\n");
    for(r = 0; r < num; r++){
        printf("%d\t%0.3e\n", r, massDelta[r]);
    }
    printf("\n");

    for(r = 0; r < num; r++){
        if(massDelta[r] < 1e-15){
            N++;
        }
    }

    printf("N = %d\n", N);
    printf("\n");

    massDelta[0] = massDelta[N + 1];
    for(l = 1; l < 5; l++){
        massDelta[l] = fabs(1 - S_finder(a, massK, massW, massD, l, N + 1, m));
        massRel[l] = massDelta[l - 1] / massDelta[l];
    }

    printf("l\tdelta\t\trel\n");
    for(l = 1; l < 5; l++){
        printf("%d\t%0.3e\t%0.3e\n", l, massDelta[l], massRel[l]);
    }
    printf("\n");

    delete [] massW;
    delete [] massD;
    delete [] massK;
    delete [] massS;
    delete [] massDelta;
    delete [] massRel;

    return 1;
}

double K_finder(double a, int r){

    double k;

    k = (r + 1) / (pow(a + 1, r + 1) - pow(a, r + 1));

    return k;
}

double S_finder(double a, double *k, double *w, double *d, int l, int r, int m){

    double s = 0;
    int i;

    for(i = 0; i < pow(2, l); i++){
        s += (1./pow(2, l)) * summF(a + i*(1./pow(2, l)), a + (i + 1)*(1./pow(2, l)), k[r], w, d, r, m);
    }

    return s;
}

double summF(double a, double b, double k, double *w, double *d, int r, int m){

    int i;
    double summ = 0;

    for(i = 0; i < m; i++){
        summ += w[i]*k*pow((a + b)/2. + d[i]*(b - a)/2., r);
    }

    return summ;
}