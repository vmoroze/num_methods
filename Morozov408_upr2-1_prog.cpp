#include <iostream>
#include <cmath>

void fill_A(double *A, int n);
void fill_B(double *b, int n, double degF);
void mult(double *A, double *c, double *res, int n);
double scal(double *a, double *b, int n);
double norm(double *a, int n);
void summ(double *a, double *b, double *res, int n);
void div(double *a, double *b, double *res, int n);
int M(double *A, double *b, double *c, double *c0, double *r, double *r0, double *p, double *p0, double *temp, int n, double eps, double *rr);

int main(void){

    double degF;
    double eps;
    int n;
    int k, i, j;
    double *A;
    double *b;
    double *c;
    double *c0;
    double *r;
    double *r0;
    double *t1d;
    int *t1t;
    double *t2d;
    int *t2t;
    double *temp;
    double *p;
    double *p0;
    double *rr;

    try{
		A = new double[65*65];
        b = new double[65];
        c = new double[65];
        c0 = new double[65];
        r = new double[65];
        r0 = new double[65];
        t1d = new double[16];
        t1t = new int[16];
        t2d = new double[16];
        t2t = new int[16];
        temp = new double[65];
        p = new double[65];
        p0 = new double[65];
        rr = new double[65];
	} catch(...){
		std::cout << "some trouble with memory" << std::endl;
		return -1;
	}

    eps = 1.000e-4;

    for(k = 3; k < 7; k++){

        n = 1;
        for(i = 0; i < k; i++){
            n *= 2;
        }

        fill_A(A, n);

        for(j = 0; j < 4; j++){

            for(i = 0; i < 65; i++){
                c0[i] = 0;
            }

            degF = j;

            fill_B(b, n, degF);

            mult(A, c0, temp, n);
            div(b, temp, r0, n);

            t1t[(k-3)*4 + j] = M(A, b, c, c0, r, r0, p, p0, temp, n, eps, rr);
            mult(A, c, temp, n);
            t1d[(k-3)*4 + j] = sqrt((1./(2*j + 1)) - 2*scal(b, c, n) + scal(temp, c, n));
        }        
    }

    eps = 1.000e-6;

    for(k = 3; k < 7; k++){

        n = 1;
        for(i = 0; i < k; i++){
            n *= 2;
        }

        fill_A(A, n);

        for(j = 0; j < 4; j++){

            for(i = 0; i < 65; i++){
                c0[i] = 0;
            }

            degF = j;

            fill_B(b, n, degF);

            mult(A, c0, temp, n);
            div(b, temp, r0, n);

            t2t[(k-3)*4 + j] = M(A, b, c, c0, r, r0, p, p0, temp, n, eps, rr);
            mult(A, c, temp, n);
            t2d[(k-3)*4 + j] = sqrt((1./(2*j + 1)) - 2*scal(b, c, n) + scal(temp, c, n));
        }        
    }

    printf("\n");
    printf("eps = %0.3e:\n", eps*100);
    printf("\n");
    printf("n/f\t1\t\t\tx\t\t\tx^2\t\t\tx^3\n");
    printf("\n");
    printf("8\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t1d[0], t1t[0], t1d[1], t1t[1], t1d[2], t1t[2], t1d[3], t1t[3]);
    printf("\n");
    printf("16\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t1d[4], t1t[4], t1d[5], t1t[5], t1d[6], t1t[6], t1d[7], t1t[7]);
    printf("\n");
    printf("32\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t1d[8], t1t[8], t1d[9], t1t[9], t1d[10], t1t[10], t1d[11], t1t[11]);
    printf("\n");
    printf("64\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t1d[12], t1t[12], t1d[13], t1t[13], t1d[14], t1t[14], t1d[15], t1t[15]);
    printf("\n");

    printf("\n");
    printf("eps = %0.3e:\n", eps);
    printf("\n");
    printf("n/f\t1\t\t\tx\t\t\tx^2\t\t\tx^3\n");
    printf("\n");
    printf("8\t%0.3e : %d\t\t%0.3e : %d\t\t\t%0.3e : %d\t\t%0.3e : %d\n", t2d[0], t2t[0], t2d[1], t2t[1], t2d[2], t2t[2], t2d[3], t2t[3]);
    printf("\n");
    printf("16\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t2d[4], t2t[4], t2d[5], t2t[5], t2d[6], t2t[6], t2d[7], t2t[7]);
    printf("\n");
    printf("32\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t2d[8], t2t[8], t2d[9], t2t[9], t2d[10], t2t[10], t2d[11], t2t[11]);
    printf("\n");
    printf("64\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t2d[12], t2t[12], t2d[13], t2t[13], t2d[14], t2t[14], t2d[15], t2t[15]);
    printf("\n");

    delete [] A;
    delete [] b;
    delete [] c;
    delete [] c0;
    delete [] r;
    delete [] r0;
    delete [] t1d;
    delete [] t1t;
    delete [] t2d;
    delete [] t2t;
    delete [] temp;
    delete [] p;
    delete [] p0;
    delete [] rr;

    return 1;
}

void fill_A(double *A, int n){

    double h = 1./n;
    int i, j;

    for(i = 0; i < n+1; i++){
        for(j = 0; j < n+1; j++){
            if(i + j == 0){
                A[(n+1)*i + j] = h/3;
                continue;
            }
            if(i + j == 2*n){
                A[(n+1)*i + j] = h/3;
                continue;
            }
            if(i == j){
                A[(n+1)*i + j] = 2*h/3;
                continue;
            }
            if(fabs(i - j) == 1){
                A[(n+1)*i + j] = h/6;
                continue;
            }
            A[(n+1)*i + j] = 0;
        }
    }

    return;
}

void fill_B(double *b, int n, double degF){

    double h = 1./n;
    int i;

    for(i = 0; i < n+1; i++){
        if(i == 0){
            b[i] = (h/6) * (pow(0, degF) + 2*pow(h/2, degF));
            continue;
        }
        if(i == n){
            b[i] = (h/6) * (pow(1, degF) + 2*pow((1 - h/2), degF));
            continue;
        }
        b[i] = (h/3) * (pow((i*h - h/2), degF) + pow(i*h, degF) + pow((i*h + h/2), degF));
    }

    return;
}

void mult(double *A, double *c, double *res, int n){

    int i, j;

    for(i = 0; i < n+1; i++){
        res[i] = 0;
    }

    for(i = 0; i < n+1; i++){
        for(j = 0; j < n+1; j++){
            res[i] += A[(n+1)*i + j] * c[j];
        }
    }

    return;
}

double scal(double *a, double *b, int n){

    int i;
    double res = 0;

    for(i = 0; i < n+1; i++){
        res += a[i] * b[i];
    }

    return res;
}

double norm(double *a, int n){

    int i;
    double res = 0;

    for(i = 0; i < n+1; i++){
        res += a[i] * a[i];
    }

    return sqrt(res);
}

void summ(double *a, double *b, double *res, int n){

    int i;

    for(i = 0; i < n+1; i++){
        res[i] = a[i] + b[i];
    }

    return;
}

void div(double *a, double *b, double *res, int n){

    int i;

    for(i = 0; i < n+1; i++){
        res[i] = a[i] - b[i];
    }

    return;
}

int M(double *A, double *b, double *c, double *c0, double *r, double *r0, double *p, double *p0, double *temp, int n, double eps, double *rr){

    int i, j;
    double aa;
    double bb;
    double norm0;

    norm0 = norm(r0, n);

    for(j = 0; j < n+1; j++){
        p0[j] = r0[j];
    }
    for(j = 0; j < n+1; j++){
        rr[j] = r0[j];
    }

    for(i = 0; norm0 / norm(rr, n) < 1./eps; i++){
        mult(A, p0, temp, n);
        aa = scal(r0, r0, n) / scal(temp, p0, n);
        for(j = 0; j < n+1; j++){
            c[j] = c0[j] + aa*p0[j];
        }
        for(j = 0; j < n+1; j++){
            r[j] = r0[j] - aa*temp[j];
        }
        bb = scal(r, r, n) / scal(r0, r0, n);
        for(j = 0; j < n+1; j++){
            p[j] = r[j] + bb*p0[j];
        }

        mult(A, c, temp, n);
        div(b, temp, rr, n);

        for(j = 0; j < n+1; j++){
            p0[j] = p[j];
            r0[j] = r[j];
            c0[j] = c[j];
        }
    }

    return i;
}