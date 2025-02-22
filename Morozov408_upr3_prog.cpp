#include <iostream>
#include <cmath>

double t0_res(double r, double s);
double sol_value(double r, double s, double *t, double *m1, double *m2);
double appr_value(double r, double s, double *t, double *WTAB, double *A1, double *A2, double *A3, int KRSUM);
int matrix_deg(double deg, double *m);
int matrix_coef(double deg, double *m, double *t, int n);
double mult(double a, double b, double c);
double fact(double d);
double mon(double c, double i, double j, double k);

int main(void){

    double *WTAB;
    double *A1;
    double *A2;
    double *A3;
    double *t;
    double *SOL;
    double *APPR;
    double *TAB;
    double *m1;
    double *m2;
    double dr, ds;
    int r, s;
    int i = 0;
    int KRSUM = 7;

    try{
		WTAB = new double[7];
        A1 = new double[7];
        A2 = new double[7];
        A3 = new double[7];
        t = new double[6];
        SOL = new double[49];
        APPR = new double[49];
        TAB = new double[49];
        m1 = new double[144];
        m2 = new double[144];
	} catch(...){
		std::cout << "some trouble with memory" << std::endl;
		return -1;
	}

    WTAB[0] = 0.225000000000000;
    WTAB[1] = WTAB[2] = WTAB[3] = 0.125939180544827;
    WTAB[4] = WTAB[5] = WTAB[6] = 0.132394152788506;
    A1[0] = A2[0] = A3[0] = 0.333333333333333;
    A1[1] = A2[2] = A3[3] = 0.797426985353087;
    A1[2] = A1[3] = A2[1] = A2[3] = A3[1] = A3[2] = 0.101286507323456;
    A1[4] = A2[5] = A3[6] = 0.059715871789770;
    A1[5] = A1[6] = A2[4] = A2[6] = A3[4] = A3[5] = 0.470142064105115;
    t[0] = 1; t[1] = 1; t[2] = 3; t[3] = 1; t[4] = 1; t[5] = 2;

    for(r = 0; r < 7; r++){
        for(s = 0; s < 7; s++){
            dr = (double)r;
            ds = (double)s;
            SOL[i] = sol_value(dr, ds, t, m1, m2);
            APPR[i] = appr_value(dr, ds, t, WTAB, A1, A2, A3, KRSUM);
            TAB[i] = fabs( (SOL[i] - APPR[i]) / SOL[i] );
            i++;
        }
    }

    printf("\n");
    printf("r/s\t0\t\t1\t\t2\t\t3\t\t4\t\t5\t\t6\n");
    printf("\n");
    for(i = 0; i < 7; i++){
        printf("%d\t%0.3e\t%0.3e\t%0.3e\t%0.3e\t%0.3e\t%0.3e\t%0.3e\n", i, TAB[i*7], TAB[i*7 + 1], TAB[i*7 + 2], TAB[i*7 + 3], TAB[i*7 + 4], TAB[i*7 + 5], TAB[i*7 + 6]);
        printf("\n");
    }

    delete [] WTAB;
    delete [] A1;
    delete [] A2;
    delete [] A3;
    delete [] t;
    delete [] SOL;
    delete [] APPR;
    delete [] TAB;
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

double sol_value(double r, double s, double *t, double *m1, double *m2){

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

double appr_value(double r, double s, double *t, double *WTAB, double *A1, double *A2, double *A3, int KRSUM){

    double appr = 0;
    int i;
    double x, y, px, py;

    for(i = 0; i < KRSUM; i++){
        x = t[0] * A1[i] + t[2] * A2[i] + t[4] * A3[i];
        y = t[1] * A1[i] + t[3] * A2[i] + t[5] * A3[i];
        px = pow(x, r);
        py = pow(y, s);
        appr += px * py * WTAB[i];
    }

    return appr;
}