#include <iostream>
#include <cmath>

double t0_res(double r, double s);
double sol_value(double r, double s, double *t, double *m1, double *m2);
double appr_value(double r, double s, double *t, double *WTAB, double *A1, double *A2, double *A3, int KRSUM, double meas);
int matrix_deg(double deg, double *m);
int matrix_coef(double deg, double *m, double *t, int n);
double mult(double a, double b, double c);
double fact(double d);
double mon(double c, double i, double j, double k);
int gen_mass(double a, double b, double c, double ***m);

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
    double ***xi;
    double ***yi;
    double *D;
    double *R;
    double *SL;
    double dr, ds;
    int r, s;
    int i, j;
    int KRSUM = 3;
    int l, N;

    try{
		WTAB = new double[3];
        A1 = new double[3];
        A2 = new double[3];
        A3 = new double[3];
        t = new double[6];
        SOL = new double[16];
        APPR = new double[16];
        TAB = new double[16];
        m1 = new double[144];
        m2 = new double[144];
        D = new double[5];
        R = new double[5];
        SL = new double[5];

        xi = new double**[257];
        for(i = 0; i < 257; i++){
            xi[i] = new double*[4];
            for(j = 0; j < 4; j++){
                xi[i][j] = new double[5];
            }
        }
        yi = new double**[257];
        for(i = 0; i < 257; i++){
            yi[i] = new double*[4];
            for(j = 0; j < 4; j++){
                yi[i][j] = new double[5];
            }
        }

	} catch(...){
		std::cout << "some trouble with memory" << std::endl;
		return -1;
	}

    WTAB[0] = WTAB[1] = WTAB[2] = 0.333333333333333;
    A1[0] = A2[1] = A3[2] = 0.000000000000000;
    A1[1] = A1[2] = A2[0] = A2[2] = A3[0] = A3[1] = 0.500000000000000;
    t[0] = 1; t[1] = 1; t[2] = 3; t[3] = 1; t[4] = 1; t[5] = 2;

    i = 0;
    for(r = 0; r < 4; r++){
        for(s = 0; s < 4; s++){
            dr = (double)r;
            ds = (double)s;
            SOL[i] = sol_value(dr, ds, t, m1, m2);
            APPR[i] = appr_value(dr, ds, t, WTAB, A1, A2, A3, KRSUM, 1);
            TAB[i] = fabs( (SOL[i] - APPR[i]) / SOL[i] );
            i++;
        }
    }

    printf("\n");
    printf("r/s\t0\t\t1\t\t2\t\t3\n");
    printf("\n");
    for(i = 0; i < 4; i++){
        printf("%d\t%0.3e\t%0.3e\t%0.3e\t%0.3e\n", i, TAB[i*4], TAB[i*4 + 1], TAB[i*4 + 2], TAB[i*4 + 3]);
        printf("\n");
    }
    printf("N = 2, r = 2, s = 1\n");
    printf("\n");
    
    r = 2;
    s = 1;
    xi[0][0][0] = xi[0][0][1] = xi[0][0][2] = xi[0][0][3] = xi[0][0][4] = 0;
    yi[0][0][0] = yi[0][0][1] = yi[0][0][2] = yi[0][0][3] = yi[0][0][4] = 0;
    xi[222][0][0] = 0;
    yi[222][0][0] = 0;
    
    gen_mass(t[0], t[2], t[4], xi);
    gen_mass(t[1], t[3], t[5], yi);

    D[0] = TAB[9];
    SL[0] = SOL[9];
    for(i = 1; i < 5; i++){
        SL[i] = 0;
    }

    for(l = 1; l < 5; l++){
        for(N = 1; N < pow(4, l) + 1; N++){
            t[0] = xi[N][1][l];
            t[1] = yi[N][1][l];
            t[2] = xi[N][2][l];
            t[3] = yi[N][2][l];
            t[4] = xi[N][3][l];
            t[5] = yi[N][3][l];
            SL[l] += appr_value((double)r, (double)s, t, WTAB, A1, A2, A3, KRSUM, 1/pow(4, l));
        }
        D[l] = fabs( (SL[0] - SL[l]) / SL[0] );
        R[l] = D[l-1] / D[l];
    }

    printf("l\t\t1\t\t2\t\t3\t\t4\n");
    printf("\n");
    printf("D[l]\t\t%0.3e\t%0.3e\t%0.3e\t%0.3e\n", D[1], D[2], D[3], D[4]);
    printf("\n");
    printf("D[l-1]/D[l]\t%0.3e\t%0.3e\t%0.3e\t%0.3e\n", R[1], R[2], R[3], R[4]);
    printf("\n");

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
    delete [] xi;
    delete [] yi;
    delete [] D;
    delete [] R;
    delete [] SL;

    for(i = 0; i < 257; i++){
        for(j = 0; j < 4; j++){
            delete [] xi[i][j];
        }
        delete [] xi[i];
    }
    delete [] xi;
    for(i = 0; i < 257; i++){
        for(j = 0; j < 4; j++){
            delete [] yi[i][j];
        }
        delete [] yi[i];
    }
    delete [] yi;

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

double appr_value(double r, double s, double *t, double *WTAB, double *A1, double *A2, double *A3, int KRSUM, double meas){

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

    return meas * appr;
}

int gen_mass(double a, double b, double c, double ***m){

    double ab, ac, bc;

    m[0][0][(int)m[222][0][0]] += 1;

    m[(int)m[0][0][(int)m[222][0][0]]][1][(int)m[222][0][0]] = a;
    m[(int)m[0][0][(int)m[222][0][0]]][2][(int)m[222][0][0]] = b;
    m[(int)m[0][0][(int)m[222][0][0]]][3][(int)m[222][0][0]] = c;

    if(m[222][0][0] == 4){
        return 1;
    }

    m[222][0][0] += 1;

    ab = (a + b) / 2;
    ac = (a + c) / 2;
    bc = (b + c) / 2;

    gen_mass(a, ab, ac, m);
    gen_mass(ab, b, bc, m);
    gen_mass(ac, bc, c, m);
    gen_mass(ac, ab, bc, m);

    m[222][0][0] -= 1;

    return 1;
}