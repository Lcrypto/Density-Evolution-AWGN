#include <stdio.h>
#include <math.h>
#include "mmatrix.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    MATRIX *P;    /*input 0*/
    MATRIX *Kin;  /*input 1*/
    MATRIX *S = matrix("S",0,0);
    MATRIX *h = matrix("h",0,0);
    MATRIX *astar = matrix("astar",0,0);
    MATRIX *Qstar = matrix("Qstar",0,0);
    MATRIX *g = matrix("g",0,0);
    MATRIX *mi = matrix("mi",1,1);
    int J,I,K;
    MWORD log2b = (1.0 / log(2));
    
    double P_ij0,P_ij1,t;
    int ii,jj,k,ap,a;
    
    
    if (nrhs != 2)
        mexErrMsgTxt("Incorrect number of input args.");
    

    P = matlab_to_matrix(prhs[0]);
    Kin = matlab_to_matrix(prhs[1]);
    J = P->nr;
    I = P->nc;
    K = (int) Kin->elt[0][0];
    
    if (K >= I) {
        Qstar = matrix("Qstar",I,I);
        for (ii=0;ii<J;ii++) {
            Qstar->elt[ii][ii] = 1;
        }
        plhs[0]= matrix_to_matlab(Qstar);
        matrix_free(g);
        matrix_free(Qstar);
        matrix_free(astar);
        matrix_free(h);
        matrix_free(S);
        matrix_free(Kin);
        matrix_free(P);
    }
    
    
    /* 
       for loops start at index 1, like Matlab.
       but subtract 1 when accessing arrays
     */
    g = matrix("g",I,I); 
    for (ap=1;ap<=I;ap++) {
        P_ij0 = 0;
        P_ij1 = 0;
        for (a=ap;a<=I;a++) {
            P_ij0 = P_ij0 + P->elt[a-1][0];
            P_ij1 = P_ij1 + P->elt[a-1][1];
            t = P_ij0 + P_ij1;
            g->elt[a-1][ap-1] = 0.5 * (P_ij0 * log2b * log(2.0 * P_ij0 / t) + P_ij1 * log2b * log(2.0 * P_ij1 / t)) + 1E-15;
        }
    }
    
    S = matrix("S",I,K);
    h = matrix("h",I,K);
    
    for (a=1;a<=1+I-K;a++) {
        S->elt[0][a-1] = g->elt[a-1][0];
    }
    
    for (k=2;k<=K;k++) {
        for (a=k;a<=k+I-K;a++) {
            S->elt[k-1][a-1] = -10000.0;
            for (ap=k-1;ap<=a-1;ap++) {
                t = S->elt[k-2][ap-1] + g->elt[a-1][ap];
                if (t > S->elt[k-1][a-1]) {
                    S->elt[k-1][a-1] = t;
                    h->elt[k-1][a-1] = ap;
                }
            }
        }
    }

    astar = matrix("astar",1,K+1);
    astar->elt[K][0] = I;
    for (k = K-1;k>=1;k--) {
        astar->elt[k][0] = h->elt[k][(int) astar->elt[k+1][0] - 1 ];
    }
    
    Qstar = matrix("Qstar",K,I);
    for (k=1;k<=K;k++) {
        for (a = astar->elt[k-1][0]+1 ; a <= astar->elt[k][0] ; a++ ) {
            Qstar->elt[a-1][k-1] = 1;
        }
    }
    
    mi->elt[0][0] = S->elt[K-1][I-1];
    plhs[0]= matrix_to_matlab(Qstar);
    plhs[1]= matrix_to_matlab(mi);

    matrix_free(mi);
    matrix_free(g);
    matrix_free(Qstar);
    matrix_free(astar);
    matrix_free(h);
    matrix_free(S);
    matrix_free(Kin);
    matrix_free(P);
    
    return;
}

