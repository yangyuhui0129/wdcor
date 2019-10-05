/*
   wdcov.c: weighted distance correlation and covariance statistics
   and wdCov test for multivariate independence

   Canhong Wen, Yuhui Yang, Quan Xiao, Meiyan Huang and Wenliang (2019)
   "Cenome-wide Association Studies of Brain Imaging Data via Weighted
    distance Correlation"
   Bioinformatics.

   Software: Maria Rizzo     mrizzo at bgsu.edu
             URL: personal.bgsu.edu/~mrizzo

*/

#include <R.h>
#include <Rmath.h>
#include <stdio.h> 


void   wdCOVtest(double *x, double *y, int *byrow, int *dims,double *gammalist,
                   double *index,  double *reps1,double *reps2,double *DCOV,
                   double *pval);
void   wdCOV(double *x, double *y, int *byrow, int *dims,double *gammalist,
               double *index, int *idx, double *DCOV);
double Akl(double **akl, double **A, int n);



/* functions in utilities.c */
extern double **alloc_matrix(int r, int c);
extern int    **alloc_int_matrix(int r, int c);
extern void   free_matrix(double **matrix, int r, int c);
extern void   free_int_matrix(int **matrix, int r, int c);
extern void   permute(int *J, int n);
extern void   roworder(double *x, int *byrow, int r, int c);
extern void   Euclidean_distance(double *x, double **Dx, int n, int d);
extern void   index_distance(double **Dx, int n, double index);
extern void   vector2matrix(double *x, double **y, int N, int d, int isroworder);
extern void   newvectorcollect(double *x, double *y, int N, int d, int r,int isroworder); 



void wdCOVtest(double *x, double *y, int *byrow, int *dims,double *gammalist,
                 double *index, double *reps1,double *reps2,
                 double *DCOV, double *pval) {
    int    i, j, k, n, n2, p, q, r, M1,M2, R,g,h,G,H,m,d,ss;
    
    int*   perm;
    double **Dx, **Dy, **A, **B;
    double dcov, V,u,dcor,dVar1,dVar2;
    double *N,*w,*t,*w1,*w2,*w3,*C,*D;
    
    n = dims[0];
    p = dims[1];
    q = dims[2];
    R = dims[3];
    ss=dims[4];
    n2 = ((double) n) * n;
    m=n*q;
    
    
    
    w = Calloc(q, double);
    w1 = Calloc(q, double);
    w2 = Calloc(q, double);
    w3 = Calloc(q, double);
    u=0.0;
    if (*byrow == FALSE) {
        /* avoid this step: use as.double(t(x)) in R */
        roworder(x, byrow, n, p);
        *byrow = FALSE;  /* false for y */
        roworder(y, byrow, n, q);
    }
    
    
    Dx = alloc_matrix(n, n);
    Euclidean_distance(x, Dx, n, p);
    index_distance(Dx, n, *index);
    A = alloc_matrix(n, n);
    Akl(Dx, A, n);
    free_matrix(Dx, n, n);
    
    for(i=0;i<q;i++){
        Dy = alloc_matrix(n, n);
        
        C = Calloc(n, double);
        newvectorcollect(y,C,n,q,i,1);
        
        Euclidean_distance(C, Dy, n, 1);
        
        index_distance(Dy, n, *index);
        
        B = alloc_matrix(n, n);
        Akl(Dy, B, n);
        
        free_matrix(Dy, n, n);
        Free(C);
        
        w[i]=0.0;
        w1[i]=0.0;
        w2[i]=0.0;
        w3[i]=0.0;
        for (k=0; k<n; k++)
            for (j=0; j<n; j++)
            {
                w1[i] += A[k][j]*B[k][j];
                w2[i] += A[k][j]*A[k][j];
                w3[i] += B[k][j]*B[k][j];
            }
        w1[i] /= n2;
        if (w1[i] > 0)
            w1[i] = sqrt(w1[i]);
        else w1[i] = 0.0;
        
        w2[i] /= n2;
        if (w2[i] > 0)
            w2[i] = sqrt(w2[i]);
        else w2[i] = 0.0;
        
        w3[i] /= n2;
        if (w3[i] > 0)
            w3[i] = sqrt(w3[i]);
        else w3[i] = 0.0;
        
        V = w2[i]*w3[i];
        if (V > DBL_EPSILON)
            w[i] = w1[i] / sqrt(V);
        else w[i] = 0.0;
        free_matrix(B, n, n);
        
    }
    
    Free(w1);Free(w2);Free(w3);
    
    N=Calloc(4,double);
    N[0]=-10000;
    N[1]=-10000;
    N[2]=-10000;
    N[3]=-10000;
    
    for(d=0;d<ss;d++){
        D = Calloc(m, double);
        t= Calloc(q, double);
        for(i=0;i<q;i++)
            t[i]=R_pow(w[i],gammalist[d]);
        
        u=0.0;
        for(i=0;i<q;i++)
            u+=R_pow(t[i],2);
        
        if(u>0)
            u=sqrt(u);
        else u=0.0;
        
        for(i=0;i<q;i++)
        {if(u>0)
            t[i]/=u;
        else t[i]=0.0;
        }
        
        for(i=0;i<q;i++)
            for(j=i;j<m;j+=q)
                D[j]=y[j]*t[i];
        Free(t);
        
        Dy = alloc_matrix(n, n);
        
        Euclidean_distance(D, Dy, n, q);
        
        Free(D);
        index_distance(Dy, n, *index);
        
        B = alloc_matrix(n, n);
        Akl(Dy, B, n);
        free_matrix(Dy, n, n);
        
        /* compute dCov(x,y), dVar(x), dVar(y) */
        for (k=0; k<4; k++)
            DCOV[k] = 0.0;
        for (k=0; k<n; k++)
            for (j=0; j<n; j++) {
                DCOV[0] += A[k][j]*B[k][j];
                DCOV[2] += A[k][j]*A[k][j];
                DCOV[3] += B[k][j]*B[k][j];
            }
        
        for (k=0; k<4; k++) {
            DCOV[k] /= n2;
            if (DCOV[k] > 0)
                DCOV[k] = sqrt(DCOV[k]);
            else DCOV[k] = 0.0;
        }
        /* compute dCor(x, y) */
        V = DCOV[2]*DCOV[3];
        if (V > DBL_EPSILON)
            DCOV[1] = DCOV[0] / sqrt(V);
        else DCOV[1] = 0.0;
        free_matrix(B, n, n);
        if(DCOV[1]>N[1])
        {N[0]=DCOV[0];N[1]=DCOV[1];N[2]=DCOV[2];N[3]=DCOV[3];}
    }
    
    DCOV[0]=N[0];DCOV[1]=N[1];DCOV[2]=N[2];DCOV[3]=N[3];
    
    Free(w);Free(N);
    
    if (R > 0) {
        /* compute the replicates */
        if (DCOV[1] > 0.0) {
            perm = Calloc(n, int);
            M1 = 0;
            M2=0;
            for (i=0; i<n; i++) perm[i] = i;
            GetRNGstate();
            for (r=0; r<R; r++) {
                permute(perm, n);
                
                w = Calloc(q, double);
                w1 = Calloc(q, double);
                w2 = Calloc(q, double);
                w3 = Calloc(q, double);
                for(i=0;i<q;i++){
                    
                    Dy = alloc_matrix(n, n);
                    
                    C = Calloc(n, double);
                    newvectorcollect(y,C,n,q,i,1);
                    
                    Euclidean_distance(C, Dy, n, 1);
                    Free(C);
                    
                    index_distance(Dy, n, *index);
                    
                    B = alloc_matrix(n, n);
                    Akl(Dy, B, n);
                    free_matrix(Dy, n, n);
                    
                    /* compute dCov(x,y), dVar(x), dVar(y) */
                    w[i]=0.0;
                    w1[i]=0.0;
                    w2[i]=0.0;
                    w3[i]=0.0;
                    for (g=0; g<n; g++)
                    {G=perm[g];
                        for (h=0; h<n; h++)
                        {
                            H=perm[h];
                            w1[i] += A[g][h]*B[G][H];
                            w2[i] += A[g][h]*A[g][h];
                            w3[i] += B[G][H]*B[G][H];
                        }
                    }
                    w1[i] /= n2;
                    if (w1[i] > 0)
                        w1[i] = sqrt(w1[i]);
                    else w1[i] = 0.0;
                    
                    w2[i] /= n2;
                    if (w2[i] > 0)
                        w2[i] = sqrt(w2[i]);
                    else w2[i] = 0.0;
                    
                    w3[i] /= n2;
                    if (w3[i] > 0)
                        w3[i] = sqrt(w3[i]);
                    else w3[i] = 0.0;
                    
                    V = w2[i]*w3[i];
                    if (V > DBL_EPSILON)
                        w[i] = w1[i] / sqrt(V);
                    else w[i] = 0.0;
                    free_matrix(B, n, n);
                    
                }
                
                Free(w1);Free(w2);Free(w3);
                
                N=Calloc(2,double);
                N[0]=-10000;
                N[1]=-10000;
                
                for(d=0;d<ss;d++){
                    D = Calloc(m, double);
                    t= Calloc(q, double);
                    for(i=0;i<q;i++)
                        t[i]=R_pow(w[i],gammalist[d]);
                    
                    u=0.0;
                    for(i=0;i<q;i++)
                        u+=R_pow(t[i],2);
                    
                    if(u>0)
                        u=sqrt(u);
                    else u=0.0;
                    
                    for(i=0;i<q;i++)
                    {if(u>0)
                        t[i]/=u;
                    else t[i]=0.0;
                    }
                    
                    for(i=0;i<q;i++)
                        for(j=i;j<m;j+=q)
                            D[j]=y[j]*t[i];
                    Free(t);
                    
                    Dy = alloc_matrix(n, n);
                    
                    Euclidean_distance(D, Dy, n, q);
                    
                    Free(D);
                    
                    index_distance(Dy, n, *index);
                    
                    B = alloc_matrix(n, n);
                    Akl(Dy, B, n);
                    free_matrix(Dy, n, n);
                    
                    dcov=0.0;
                    dVar1=0.0;
                    dVar2=0.0;
                    dcor=0.0;
                    for (g=0; g<n; g++)
                    {G=perm[g];
                        for (h=0; h<n; h++) {
                            H=perm[h];
                            dcov += A[g][h]*B[G][H];
                            dVar1 +=A[g][h]*A[g][h];
                            dVar2 +=B[G][H]*B[G][H];
                        }
                    }
                    
                    dcov /= n2;
                    if (dcov > 0)
                        dcov = sqrt(dcov);
                    else dcov = 0.0;
                    dVar1 /= n2;
                    if (dVar1 > 0)
                        dVar1 = sqrt(dVar1);
                    else dVar1 = 0.0;
                    
                    dVar2 /= n2;
                    if (dVar2 > 0)
                        dVar2 = sqrt(dVar2);
                    else dVar2 = 0.0;
                    
                    V = dVar1*dVar2;
                    if (V > DBL_EPSILON)
                        dcor =dcov / sqrt(V);
                    else dcor = 0.0;
                    
                    free_matrix(B, n, n);
                    
                    if(dcor>N[1])
                    {N[0]=dcov;N[1]= dcor;}
                }
                dcov=N[0];dcor=N[1];
                Free(w);Free(N);
                reps1[r] = dcov;
                reps2[r] = dcor;
                if (dcov >= DCOV[1]) M1++;
                if (dcor >= DCOV[1]) M2++;
            }
            pval[0] = (double) (M1+1) / (double) (R+1);
            pval[1] = (double) (M2+1) / (double) (R+1);
            PutRNGstate();
            Free(perm);
        }
        else {
            pval[0] = 1.0;
            pval[1] = 1.0;
        }
    }
    
    free_matrix(A, n, n);
    return;
}


void wdCOV(double *x, double *y, int *byrow, int *dims,double *gammalist,
             double *index, int *idx, double *DCOV) {
    int i,j,n,p,q,k,n2,m,d,ss;
    double **Dx, **Dy,**A,**B;
    double u,V;
    double *w,*t,*w1,*w2,*w3,*C,*D,*N;
    
    n = dims[0];
    p = dims[1];//1
    q = dims[2];
    ss=dims[3];
    n2=((double) n) * n;
    m=n*q;
    
    w = Calloc(q, double);
    w1 = Calloc(q, double);
    w2 = Calloc(q, double);
    w3 = Calloc(q, double);
    u=0.0;
    
    
    if (*byrow == FALSE) {
        /* avoid this step: use as.double(t(x)) in R */
        roworder(x, byrow, n, p);
        *byrow = FALSE;  /* false for y */
        roworder(y, byrow, n, q);
    }
    
    Dx = alloc_matrix(n, n);
    Euclidean_distance(x, Dx, n, p);
    index_distance(Dx, n, *index);
    A = alloc_matrix(n, n);
    Akl(Dx, A, n);
    free_matrix(Dx, n, n);
    
    for(i=0;i<q;i++){
        Dy = alloc_matrix(n, n);
        
        C = Calloc(n, double);
        newvectorcollect(y,C,n,q,i,1);
        
        Euclidean_distance(C, Dy, n, 1);
        
        index_distance(Dy, n, *index);
        
        B = alloc_matrix(n, n);
        
        Akl(Dy, B, n);
        
        free_matrix(Dy, n, n);
        Free(C);
        
        
        /* compute dCov(x,y), dVar(x), dVar(y) */
        w[i]=0.0;
        w1[i]=0.0;
        w2[i]=0.0;
        w3[i]=0.0;
        for (k=0; k<n; k++)
            for (j=0; j<n; j++)
            {
                w1[i] += A[k][j]*B[k][j];
                w2[i] += A[k][j]*A[k][j];
                w3[i] += B[k][j]*B[k][j];
            }
        w1[i] /= n2;
        if (w1[i] > 0)
            w1[i] = sqrt(w1[i]);
        else w1[i] = 0.0;
        
        w2[i] /= n2;
        if (w2[i] > 0)
            w2[i] = sqrt(w2[i]);
        else w2[i] = 0.0;
        
        w3[i] /= n2;
        if (w3[i] > 0)
            w3[i] = sqrt(w3[i]);
        else w3[i] = 0.0;
        
        V = w2[i]*w3[i];
        if (V > DBL_EPSILON)
            w[i] = w1[i] / sqrt(V);
        else w[i] = 0.0;
        
        free_matrix(B, n, n);
        
    }
    
    Free(w1);Free(w2);Free(w3);
    
    N=Calloc(4,double);
    N[0]=-10000;
    N[1]=-10000;
    N[2]=-10000;
    N[3]=-10000;
    
    for(d=0;d<=ss;d++){
        D = Calloc(m, double);
        t= Calloc(q, double);
        for(i=0;i<q;i++)
            t[i]=R_pow(w[i],gammalist[d]);
        
        u=0.0;
        for(i=0;i<q;i++)
            u+=R_pow(t[i],2);
        
        if(u>0)
            u=sqrt(u);
        else u=0.0;
        
        for(i=0;i<q;i++)
        {if(u>0)
            t[i]/=u;
        else t[i]=0.0;
        }
        
        for(i=0;i<q;i++)
            for(j=i;j<m;j+=q)
                D[j]=y[j]*t[i];
        Free(t);
        
        Dy = alloc_matrix(n, n);
        
        Euclidean_distance(D, Dy, n, q);
        
        Free(D);
        
        index_distance(Dy, n, *index);
        
        B = alloc_matrix(n, n);
        Akl(Dy, B, n);
        free_matrix(Dy, n, n);
        
        /* compute dCov(x,y), dVar(x), dVar(y) */
        for (k=0; k<4; k++)
            DCOV[k] = 0.0;
        for (k=0; k<n; k++)
            for (j=0; j<n; j++) {
                DCOV[0] += A[k][j]*B[k][j];
                DCOV[2] += A[k][j]*A[k][j];
                DCOV[3] += B[k][j]*B[k][j];
            }
        
        for (k=0; k<4; k++) {
            DCOV[k] /= n2;
            if (DCOV[k] > 0)
                DCOV[k] = sqrt(DCOV[k]);
            else DCOV[k] = 0.0;
        }
        /* compute dCor(x, y) */
        V = DCOV[2]*DCOV[3];
        if (V > DBL_EPSILON)
            DCOV[1] = DCOV[0] / sqrt(V);
        else DCOV[1] = 0.0;
        
        free_matrix(B, n, n);
        
        if(DCOV[1]>N[1])
        {N[0]=DCOV[0];N[1]=DCOV[1];N[2]=DCOV[2];N[3]=DCOV[3];}
    }
    
    free_matrix(A, n, n);
    DCOV[0]=N[0];DCOV[1]=N[1];DCOV[2]=N[2];DCOV[3]=N[3];
    
    Free(w); Free(N);
    
    return;
}


double Akl(double **akl, double **A, int n) {
    /* -computes the A_{kl} or B_{kl} distances from the
        distance matrix (a_{kl}) or (b_{kl}) for dCov, dCor, dVar
        dCov = mean(Akl*Bkl), dVar(X) = mean(Akl^2), etc.
    */
    int j, k;
    double *akbar;
    double abar;

    akbar = Calloc(n, double);
    abar = 0.0;
    for (k=0; k<n; k++) {
        akbar[k] = 0.0;
        for (j=0; j<n; j++) {
            akbar[k] += akl[k][j];
        }
        abar += akbar[k];
        akbar[k] /= (double) n;
    }
    abar /= (double) (n*n);

    for (k=0; k<n; k++)
        for (j=k; j<n; j++) {
            A[k][j] = akl[k][j] - akbar[k] - akbar[j] + abar;
            A[j][k] = A[k][j];
        }
    Free(akbar);
    return(abar);
}



