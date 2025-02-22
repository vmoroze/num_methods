#include <iostream>
#include <cmath>

double t0_res(double r, double s);
double value(double r, double s, double *t, double *m1, double *m2);
int matrix_deg(double deg, double *m);
int matrix_coef(double deg, double *m, double *t, int n);
double mult(double a, double b, double c);
double fact(double d);
double mon(double c, double i, double j, double k);

int main(void){

    double *rs;
    double *t0;
    double *t3;
    double *otvet1;
    double *otvet2;
    double *otvet3;
    double *otvet4;
    double *m1;
    double *m2;
    int i;

    try{
		rs = new double[6];
        t0 = new double[6];
        t3 = new double[6];
        otvet1 = new double[3];
        otvet2 = new double[3];
        otvet3 = new double[3];
        otvet4 = new double[3];
        m1 = new double[144];
        m2 = new double[144];
	} catch(...){
		std::cout << "some trouble with memory" << std::endl;
		return -1;
	}

    rs[0] = 0; rs[1] = 3; rs[2] = 5; rs[3] = 0; rs[4] = 3; rs[5] = 4;
    t0[0] = 0; t0[1] = 0; t0[2] = 0; t0[3] = 1; t0[4] = 2; t0[5] = 0;
    t3[0] = 1; t3[1] = 1; t3[2] = 3; t3[3] = 1; t3[4] = 1; t3[5] = 2;

    for(i = 0; i < 144; i++){
        m1[i] = 0;
        m2[i] = 0;
    }

    for(i = 0; i < 3; i++){
        otvet1[i] = t0_res(rs[2*i], rs[2*i+1]);
    }
    for(i = 0; i < 3; i++){
        otvet2[i] = value(rs[2*i], rs[2*i+1], t0, m1, m2);
    }
    for(i = 0; i < 3; i++){
        otvet3[i] = fabs(otvet1[i] - otvet2[i]);
    }
    for(i = 0; i < 3; i++){
        otvet4[i] = value(rs[2*i], rs[2*i+1], t3, m1, m2);
    }

    printf("\n");
    printf("simple\t\tcomplex\t\tdifference\tanswer\n");
    for(i = 0; i < 3; i++){
        printf("%0.3e\t%0.3e\t%0.3e\t%0.3e\n", otvet1[i], otvet2[i], otvet3[i], otvet4[i]);
    }
    printf("\n");

    delete [] rs;
    delete [] t0;
    delete [] t3;
    delete [] otvet1;
    delete [] otvet2;
    delete [] otvet3;
    delete [] otvet4;
    delete [] m1;
    delete [] m2;

    return 1;
}

double t0_res(double r, double s){

    double res = 0;
    int k;

    for(k = 0; k < s+2; k++){
        res += (fact(s+1) / (fact(s+1-k) * fact(k))) * pow(-1, k) * (1 / (r+k+1)); 
    }
    res *= ( pow(2, r+1) / (s+1) );

    return res;
}

double value(double r, double s, double *t, double *m1, double *m2){

    double res = 0;
    int i, j;
    double detB;

    matrix_deg(r, m1);
    matrix_deg(s, m2);
    matrix_coef(r, m1, t, 0);
    matrix_coef(s, m2, t, 1);

    detB = fabs( t[0]*t[3] + t[4]*t[1] + t[2]*t[5] - t[3]*t[4] - t[2]*t[1] - t[0]*t[5] );

    for(i = 0; i < (r + 2.)*(r + 1.)/2; i++){
        for(j = 0; j < (s + 2.)*(s + 1.)/2; j++){
            res += mon(m1[4*i + 3] * m2[4*j + 3],
                        m1[4*i] + m2[4*j],
                        m1[4*i + 1] + m2[4*j + 1],
                        m1[4*i + 2] + m2[4*j + 2]);
        }
    }

    res *= detB;

    return res;
}

int matrix_deg(double deg, double *m){

    int i;

    m[0] = deg; m[1] = 0; m[2] = 0;
    for(i = 1; (m[4*(i-1)] != 0) || (m[4*(i-1) + 1] != 0); i++){
        if(m[4*(i-1) + 1] == 0){
            m[4*i] = m[4*(i-1)] - 1;
            m[4*i + 1] = deg - m[4*i];
            m[4*i + 2] = 0;
        } else{
            m[4*i + 1] = m[4*(i-1) + 1] - 1;
            m[4*i + 2] = m[4*(i-1) + 2] + 1;
            m[4*i] = m[4*(i-1)];
        }
    }

    return 1;
}

int matrix_coef(double deg, double *m, double *t, int n){

    int i;

    for(i = 0; i < (deg + 2.)*(deg + 1.)/2; i++){
        m[4*i + 3] = mult(m[4*i], m[4*i + 1], m[4*i + 2])*
        pow(t[n], m[4*i]) * pow(t[2 + n], m[4*i + 1]) * pow(t[4 + n], m[4*i + 2]);
        //printf("%lf %lf %lf %lf\n", m[4*i], m[4*i + 1], m[4*i + 2], m[4*i + 3]);
    }
    //printf("\n");
    
    return 1;
}

double mult(double a, double b, double c){

    double res;

    res = fact(a + b + c) / (fact(a) * fact(b) * fact(c));

    return res;
}

double fact(double d){

    double res = 1;
    int i;

    for(i = 2; i < d + 1; i++){
        res *= i;
    }

    return res;
}

double mon(double c, double i, double j, double k){

    double res;

    res = c * fact(i) * fact(j) * fact(k) / fact(i + j + k + 2);

    return res;
}