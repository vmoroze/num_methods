#include <iostream>
#include <cmath>

void fill_A(double *A, int *II, int *JJ, int n);
void fill_B(double *b, int n, double degF, double d1, double d2, double d3);
void mult(double *A, int *II, int *JJ, double *c, double *res, int n);
double fun(double x, double y, double degF, double d1, double d2, double d3);
double scal(double *a, double *b, int n);
double norm(double *a, int n);
void summ(double *a, double *b, double *res, int n);
void div(double *a, double *b, double *res, int n);
int M(double *A, int *II, int *JJ, double *b, double *c, double *c0, double *r, double *r0, double *p, double *p0, double *temp, int n, double eps, double *rr);

int main(void){

    double degF;
    double eps;
    int n;
    int k, i, j;
    double d1, d2, d3;
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
    int *II;
    int *JJ;

    try{
        II = new int[43489];
        JJ = new int[43489]; 
		A = new double[43489];
        b = new double[6305];
        c = new double[6305];
        c0 = new double[6305]; 
        r = new double[6305];
        r0 = new double[6305];
        t1d = new double[12];
        t1t = new int[12];
        t2d = new double[12];
        t2t = new int[12];
        temp = new double[6305]; 
        p = new double[6305];
        p0 = new double[6305];
        rr = new double[6305];
	} catch(...){
		std::cout << "some trouble with memory" << std::endl;
		return -1;
	}

    d1 = d2 = d3 = 1.;

    eps = 1.000e-4;

    for(k = 3; k < 7; k++){

        n = 1;
        for(i = 0; i < k; i++){
            n *= 2;
        }

        fill_A(A, II, JJ, n);

        for(j = 0; j < 3; j++){

            for(i = 0; i < 6305; i++){
                c0[i] = 0;
            }

            degF = j;

            fill_B(b, n, degF, d1, d2, d3);

            mult(A, II, JJ, c0, temp, n);
            div(b, temp, r0, n);

            t1t[(k-3)*3 + j] = M(A, II, JJ, b, c, c0, r, r0, p, p0, temp, n, eps, rr);
            mult(A, II, JJ, c, temp, n);
            if(j == 0){
                t1d[(k-3)*3 + j] = sqrt((3*d1*d1)/2. - 2*scal(b, c, n) + scal(temp, c, n));
            }
            if(j == 1){
                t1d[(k-3)*3 + j] = sqrt((15*d1*d1 + 17*d1*d2 + 7*d2*d2)/12. - 2*scal(b, c, n) + scal(temp, c, n));
            }
            if(j == 2){
                t1d[(k-3)*3 + j] = sqrt((126*d1*d1 + 74*d1*d2 + 22*d2*d2 + 129*d1*d3 + 49*d2*d3 + 37*d3*d3)/60. - 2*scal(b, c, n) + scal(temp, c, n));
            }
        }
    }

    eps = 1.000e-6;

    for(k = 3; k < 7; k++){

        n = 1;
        for(i = 0; i < k; i++){
            n *= 2;
        }

        fill_A(A, II, JJ, n);

        for(j = 0; j < 3; j++){

            for(i = 0; i < 6305; i++){
                c0[i] = 0;
            }

            degF = j;

            fill_B(b, n, degF, d1, d2, d3);

            mult(A, II, JJ, c0, temp, n);
            div(b, temp, r0, n);

            t2t[(k-3)*3 + j] = M(A, II, JJ, b, c, c0, r, r0, p, p0, temp, n, eps, rr);
            mult(A, II, JJ, c, temp, n);
            if(j == 0){
                t2d[(k-3)*3 + j] = sqrt((3*d1*d1)/2. - 2*scal(b, c, n) + scal(temp, c, n));
            }
            if(j == 1){
                t2d[(k-3)*3 + j] = sqrt((15*d1*d1 + 17*d1*d2 + 7*d2*d2)/12. - 2*scal(b, c, n) + scal(temp, c, n));
            }
            if(j == 2){
                t2d[(k-3)*3 + j] = sqrt((126*d1*d1 + 74*d1*d2 + 22*d2*d2 + 129*d1*d3 + 49*d2*d3 + 37*d3*d3)/60. - 2*scal(b, c, n) + scal(temp, c, n));
            }
        }
    }

    printf("\n");
    printf("eps = %0.3e:\n", eps*100);
    printf("\n");
    printf("n/deg\t1\t\t\t2\t\t\t3\n");
    printf("\n");
    printf("8\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t1d[0], t1t[0], t1d[1], t1t[1], t1d[2], t1t[2]);
    printf("\n");
    printf("16\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t1d[3], t1t[3], t1d[4], t1t[4], t1d[5], t1t[5]);
    printf("\n");
    printf("32\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t1d[6], t1t[6], t1d[7], t1t[7], t1d[8], t1t[8]);
    printf("\n");
    printf("64\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t1d[9], t1t[9], t1d[10], t1t[10], t1d[11], t1t[11]);
    printf("\n");

    printf("\n");
    printf("eps = %0.3e:\n", eps);
    printf("\n");
    printf("n/deg\t1\t\t\t2\t\t\t3\n");
    printf("\n");
    printf("8\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t2d[0], t2t[0], t2d[1], t2t[1], t2d[2], t2t[2]);
    printf("\n");
    printf("16\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t2d[3], t2t[3], t2d[4], t2t[4], t2d[5], t2t[5]);
    printf("\n");
    printf("32\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t2d[6], t2t[6], t2d[7], t2t[7], t2d[8], t2t[8]);
    printf("\n");
    printf("64\t%0.3e : %d\t\t%0.3e : %d\t\t%0.3e : %d\n", t2d[9], t2t[9], t2d[10], t2t[10], t2d[11], t2t[11]);
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
    delete [] II;
    delete [] JJ;

    return 1;
}

void fill_A(double *A, int *II, int *JJ, int n){

    double h = 1./n; 
    int i, j;
    int k = 0;
    int ots = 0;
    int ots1 = 0;

    for(i = 0; i < n+1; i++){
        for(j = 0; j < n+i+1; j++){

            if(i == 0){
                if(j == 0 || j == n){
                    II[k] = JJ[k] = j;
                    A[k] = 4*((h*h)/24.);
                    k++;
                } else{
                    II[k] = JJ[k] = j;
                    A[k] = 6*((h*h)/24.);
                    k++;
                }
                if(j != n){
                    II[k] = j;
                    JJ[k] = j + 1;
                    A[k] = 1*((h*h)/24.);
                    k++;
                    II[k] = j + 1;
                    JJ[k] = j;
                    A[k] = 1*((h*h)/24.);
                    k++;
                }
                continue;
            }

            if(i == n){
                if(j == 0 || j == n + ots1){
                    II[k] = JJ[k] = ots + j;
                    A[k] = 2*((h*h)/24.);
                    k++;
                    if(j == 0){
                        II[k] = j + ots - n - i;
                    } else{
                        II[k] = j + ots - n - i - 1;
                    }
                    JJ[k] = j + ots;
                    A[k] = 1*((h*h)/24.);
                    k++;
                    if(j == 0){
                        JJ[k] = j + ots - n - i;
                    } else{
                        JJ[k] = j + ots - n - i - 1;
                    }
                    II[k] = j + ots;
                    A[k] = 1*((h*h)/24.);
                    k++;
                } else{
                    II[k] = JJ[k] = ots + j;
                    A[k] = 6*((h*h)/24.);
                    k++;
                    II[k] = j + ots - n - i;
                    JJ[k] = j + ots;
                    A[k] = 2*((h*h)/24.);
                    k++;
                    II[k] = j + ots - n - i - 1;
                    JJ[k] = j + ots;
                    A[k] = 2*((h*h)/24.);
                    k++;
                    JJ[k] = j + ots - n - i;
                    II[k] = j + ots;
                    A[k] = 2*((h*h)/24.);
                    k++;
                    JJ[k] = j + ots - n - i - 1;
                    II[k] = j + ots;
                    A[k] = 2*((h*h)/24.);
                    k++;
                }
                if(j != n+i){
                    II[k] = ots + j;
                    JJ[k] = ots + j + 1;
                    A[k] = 1*((h*h)/24.);
                    k++;
                    II[k] = ots + j + 1;
                    JJ[k] = ots + j;
                    A[k] = 1*((h*h)/24.);
                    k++;
                }
                continue;
            }

            if(j == 0 || j == n + ots1){
                II[k] = JJ[k] = ots + j;
                A[k] = 6*((h*h)/24.);
                k++;
                if(j == 0){
                    II[k] = j + ots - n - i;
                } else{
                    II[k] = j + ots - n - i - 1;
                }
                JJ[k] = j + ots;
                A[k] = 1*((h*h)/24.);
                k++;
                if(j == 0){
                    JJ[k] = j + ots - n - i;
                } else{
                    JJ[k] = j + ots - n - i - 1;
                }
                II[k] = j + ots;
                A[k] = 1*((h*h)/24.);
                k++;
            } else{
                II[k] = JJ[k] = ots + j;
                A[k] = 12*((h*h)/24.);
                k++;
                II[k] = j + ots - n - i;
                JJ[k] = j + ots;
                A[k] = 2*((h*h)/24.);
                k++;
                II[k] = j + ots - n - i - 1;
                JJ[k] = j + ots;
                A[k] = 2*((h*h)/24.);
                k++;
                JJ[k] = j + ots - n - i;
                II[k] = j + ots;
                A[k] = 2*((h*h)/24.);
                k++;
                JJ[k] = j + ots - n - i - 1;
                II[k] = j + ots;
                A[k] = 2*((h*h)/24.);
                k++;
            }
            if(j != n+i){
                II[k] = ots + j;
                JJ[k] = ots + j + 1;
                A[k] = 2*((h*h)/24.);
                k++;
                II[k] = ots + j + 1;
                JJ[k] = ots + j;
                A[k] = 2*((h*h)/24.);
                k++;
            }
        }
        ots = ots + 1 + i + n;
        ots1++;
    }

    return;
}

void fill_B(double *b, int n, double degF, double d1, double d2, double d3){

    double h = 1./n;
    int xi, yi;

    for(yi = 0; yi < n+1; yi++){
        for(xi = 0; xi < n+1+yi; xi++){
            if(xi == 0 && yi == 0){
                b[ (n + 1)*yi + xi + 1 + ((yi + 1)*(yi - 2))/2 ] = 
                2*((h*h)/24.) * (fun(h/2, 0, degF, d1, d2, d3) + 2*fun(h/2, h/2, degF, d1, d2, d3) + fun(0, h/2, degF, d1, d2, d3));
                continue;
            }
            if(xi == n && yi == 0){
                b[ (n + 1)*yi + xi + 1 + ((yi + 1)*(yi - 2))/2 ] = 
                2*((h*h)/24.) * (fun(1+h/2, h/2, degF, d1, d2, d3) + 2*fun(1, h/2, degF, d1, d2, d3) + fun(1-h/2, 0, degF, d1, d2, d3));
                continue;
            }
            if(xi == 0 && yi == n){
                b[ (n + 1)*yi + xi + 1 + ((yi + 1)*(yi - 2))/2 ] = 
                2*((h*h)/24.) * (fun(h/2, 1, degF, d1, d2, d3) + fun(0, 1-h/2, degF, d1, d2, d3));
                continue;
            }
            if(xi == 2*n && yi == n){
                b[ (n + 1)*yi + xi + 1 + ((yi + 1)*(yi - 2))/2 ] = 
                2*((h*h)/24.) * (fun(2-h/2, 1, degF, d1, d2, d3) + fun(2-h/2, 1-h/2, degF, d1, d2, d3));
                continue;
            }
            if(yi == 0){
                b[ (n + 1)*yi + xi + 1 + ((yi + 1)*(yi - 2))/2 ] = 
                2*((h*h)/24.) * (fun(h*xi+h/2, 0, degF, d1, d2, d3) + 2*fun(h*xi+h/2, h/2, degF, d1, d2, d3) + 
                2*fun(h*xi, h/2, degF, d1, d2, d3) + fun(h*xi-h/2, 0, degF, d1, d2, d3));
                continue;
            }
            if(xi == 0){
                b[ (n + 1)*yi + xi + 1 + ((yi + 1)*(yi - 2))/2 ] = 
                2*((h*h)/24.) * (fun(0, h*yi-h/2, degF, d1, d2, d3) + 2*fun(h/2, h*yi, degF, d1, d2, d3) + 
                2*fun(h/2, h*yi+h/2, degF, d1, d2, d3) + fun(0, h*yi+h/2, degF, d1, d2, d3));
                continue;
            }
            if(yi == n){
                b[ (n + 1)*yi + xi + 1 + ((yi + 1)*(yi - 2))/2 ] = 
                2*((h*h)/24.) * (fun(h*xi-h/2, 1, degF, d1, d2, d3) + 2*fun(h*xi-h/2, 1-h/2, degF, d1, d2, d3) + 
                2*fun(h*xi, 1-h/2, degF, d1, d2, d3) + fun(h*xi+h/2, 1, degF, d1, d2, d3));
                continue;
            }
            if(xi == yi + n){
                b[ (n + 1)*yi + xi + 1 + ((yi + 1)*(yi - 2))/2 ] = 
                2*((h*h)/24.) * (fun(h*xi+h/2, h*yi+h/2, degF, d1, d2, d3) + 2*fun(h*xi, h*yi+h/2, degF, d1, d2, d3) + 
                2*fun(h*xi-h/2, h*yi, degF, d1, d2, d3) + fun(h*xi-h/2, h*yi-h/2, degF, d1, d2, d3));
                continue;
            }
            b[ (n + 1)*yi + xi + 1 + ((yi + 1)*(yi - 2))/2 ] = 
            4*((h*h)/24.) * (fun(h*xi+h/2, h*yi, degF, d1, d2, d3) + fun(h*xi+h/2, h*yi+h/2, degF, d1, d2, d3) + fun(h*xi, h*yi+h/2, degF, d1, d2, d3) + 
            fun(h*xi-h/2, h*yi, degF, d1, d2, d3) + fun(h*xi-h/2, h*yi-h/2, degF, d1, d2, d3) + fun(h*xi, h*yi-h/2, degF, d1, d2, d3));
        }
    }

    return;
}

void mult(double *A, int *II, int *JJ, double *c, double *res, int n){

    int i;

    for(i = 0; i < 1.5*n*n + 2.5*n + 1; i++){
        res[i] = 0;
    }

    for(i = 0; i < 10.5*n*n + 7.5*n + 1; i++){
        res[II[i]] += A[i] * c[JJ[i]];
    }

    return;
}

double scal(double *a, double *b, int n){

    int i;
    double res = 0;

    for(i = 0; i < 1.5*n*n + 2.5*n + 1; i++){
        res += a[i] * b[i];
    }

    return res;
}

double norm(double *a, int n){

    int i;
    double res = 0;

    for(i = 0; i < 1.5*n*n + 2.5*n + 1; i++){
        res += a[i] * a[i];
    }

    return sqrt(res);
}

void summ(double *a, double *b, double *res, int n){

    int i;

    for(i = 0; i < 1.5*n*n + 2.5*n + 1; i++){
        res[i] = a[i] + b[i];
    }

    return;
}

void div(double *a, double *b, double *res, int n){

    int i;

    for(i = 0; i < 1.5*n*n + 2.5*n + 1; i++){
        res[i] = a[i] - b[i];
    }

    return;
}

int M(double *A, int *II, int *JJ, double *b, double *c, double *c0, double *r, double *r0, double *p, double *p0, double *temp, int n, double eps, double *rr){

    int i, j;
    double aa;
    double bb;
    double norm0;

    norm0 = norm(r0, n);

    for(j = 0; j < 1.5*n*n + 2.5*n + 1; j++){
        p0[j] = r0[j];
    }
    for(j = 0; j < 1.5*n*n + 2.5*n + 1; j++){
        rr[j] = r0[j];
    }

    for(i = 0; norm0 / norm(rr, n) < 1./eps; i++){
        mult(A, II, JJ, p0, temp, n);
        aa = scal(r0, r0, n) / scal(temp, p0, n);
        for(j = 0; j < 1.5*n*n + 2.5*n + 1; j++){
            c[j] = c0[j] + aa*p0[j];
        }
        for(j = 0; j < 1.5*n*n + 2.5*n + 1; j++){
            r[j] = r0[j] - aa*temp[j];
        }
        bb = scal(r, r, n) / scal(r0, r0, n);
        for(j = 0; j < 1.5*n*n + 2.5*n + 1; j++){
            p[j] = r[j] + bb*p0[j];
        }

        mult(A, II, JJ, c, temp, n);
        div(b, temp, rr, n);

        for(j = 0; j < 1.5*n*n + 2.5*n + 1; j++){
            p0[j] = p[j];
            r0[j] = r[j];
            c0[j] = c[j];
        }
    }

    return i;
}

double fun(double x, double y, double degF, double d1, double d2, double d3){

    if(degF == 0){
        return d1;
    }
    if(degF == 1){
        return d1*x + d2*y;
    }

    return d1*x*x + d2*y*y + d3*x*y;
}