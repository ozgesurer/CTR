#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/time.h>
#include <math.h>
#include <stdbool.h>
#include "R.h"
#include "Rmath.h"
#include "math.h"
#include <R_ext/Utils.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <sys/time.h>
#include <time.h>

double *base_arr;
static int comparCTR (const void *a, const void *b){
  int aa = *((int *) a), bb = *((int *) b);
  if (base_arr[aa] < base_arr[bb])
    return -1;
  else if (base_arr[aa] == base_arr[bb])
    return 0;
  else
    return 1;
}
void reverseArrayCTR(int arr[], int start, int end){
  int temp;
  while (start < end)
  {
    temp = arr[start];
    arr[start] = arr[end];
    arr[end] = temp;
    start++;
    end--;
  }
}

void CTRwithoutCV(double *X, double *Y, int *dim, int *obs, int *K, double *Z,
                  double *Rsq, double *betaTree, double *isSplit, double *GroupTree, double *JTree, int *nongr){

  int i, j, k, l, iter, iterNext, index1, index2, t_j;
  int ngroups = *K; int n = *obs; int p = *dim; int np = p*n; int one = 1;
  double d1 = 1.0; double d0 = 0.0; double dm1 = -1.0;
  int ngr = ngroups + 1; int nngr = n*ngr;
  const double *cd1 = &d1; const double *cd0 = &d0; const double *cdm1 = &dm1; const int *ci1 = &one; const int *cngr = &ngr; const int *cngroups = &ngroups; const int *cn = &n;

  double *Yhat = (double *) malloc(n*sizeof(double));
  double *E = (double *) malloc(np*sizeof(double));
  double *w = (double *) malloc(p*sizeof(double)); double *u = (double *) malloc(p*sizeof(double));
  double *r = (double *) malloc(p*sizeof(double)); double *v = (double *) malloc(p*sizeof(double));
  int *Groups = (int *) malloc(p*sizeof(int)); int *J = (int *) malloc(p*sizeof(int));
  double *z = (double *) malloc(n*sizeof(double)); double *e_z = (double *) malloc(n*sizeof(double));
  double *A = (double *) malloc(p*ngroups*sizeof(double));
  double *abest1 = (double *) malloc(p*sizeof(double)); double *abest2 = (double *) malloc(p*sizeof(double));
  double *aNew = (double *) malloc(p*sizeof(double)); double *zbest = (double *) malloc(n*sizeof(double));
  double *XY = (double *) malloc(p*sizeof(double)); double *ZY = (double *) malloc(ngr*sizeof(double));
  double *ZOld = (double *) malloc(nngr*sizeof(double));
  double *ZZinv = (double *) malloc(ngr*ngr*sizeof(double)); double *Zz = (double *) malloc(ngr*sizeof(double));
  double *ZZinvZz = (double *) malloc(ngr*sizeof(double)); double *ZZinvZY = (double *) malloc(ngr*sizeof(double));
  double *alpha = (double *) malloc(ngroups*sizeof(double)); double *beta = (double *) malloc(p*sizeof(double));
  double *colsum = (double *) malloc(p*sizeof(double));
  double ZzZZinvZz, denum, invvalue, zz, Dbestinv, abestXY, SSE, SST;
  double xvalue, yvalue, uvalue, wvalue, evalue, vvalue, zvalue, xyvalue, qvalue;
  double sumY = 0; int groupsize = 0;
  double Rbest, Nbest, Dbest, R, Dold, Dnew, Nold, Nnew, Nnewsq;
  int splitGroup, splitPoint, splitsize;

  memset(ZOld, 0, nngr*sizeof(double)); memset(ZZinv, 0, ngr*ngr*sizeof(double));
  memset(Zz, 0, ngr*sizeof(double)); memset(ZZinvZz, 0, ngr*sizeof(double));
  memset(ZY, 0, ngr*sizeof(double)); memset(e_z, 0, n*sizeof(double));
  memset(A, 0, p*ngroups*sizeof(double)); memset(v, 0, p*sizeof(double));

  for (j = 0; j < n; j++){
    Z[j] = 1;
    sumY = sumY + Y[j];
  }
  ZY[0] = sumY; ZZinv[0] = (double) 1/n;
  sumY = (double) sumY/n;
  SST = 0;
  for (j = 0; j < n; j++)
    SST = SST + pow(Y[j] - sumY, 2);

  for (i = 0; i < p; i++){
    colsum[i] = 0;
    for (j = 0; j < n; j++){
      colsum[i] = colsum[i] + X[j + n*i];
    }
    colsum[i] = colsum[i]/n;
  }

  for (i = 0; i < p; i++){
    index1 = i*n;
    uvalue = 0; wvalue = 0; xyvalue = 0;
    XY[i] = 0; Groups[i] = 0; J[i] = i;
    for(j = 0; j < n; j++){
      xvalue = X[index1 + j];
      evalue = xvalue - colsum[i];
      E[index1 + j] = evalue;
      yvalue = Y[j];
      xyvalue =  xyvalue + xvalue*yvalue;
      uvalue = uvalue + evalue*evalue;
      wvalue = wvalue + evalue*yvalue;
    }
    u[i] = uvalue; w[i] = wvalue; XY[i] = xyvalue;
    if(u[i] == 0)
      r[i] = 0;
    else{
      r[i] = copysign(1.0, w[i])*pow(w[i], 2)/u[i];
    }
  }

  for (iter = 0; iter < ngroups; iter++){

    iterNext = iter + 1;
    Rbest = 0; Nbest = 0; Dbest = 0; splitGroup = 0; splitPoint = 0; splitsize = 0;
    memset(abest1, 0, p*sizeof(double)); memset(abest2, 0, p*sizeof(double)); memset(zbest, 0, n*sizeof(double));

    for (k = 0; k < iterNext; k++){
      int ids[p]; groupsize = 0; int startGr = 0;
      for(i = 0; i < p; i++){
        if(Groups[i] == k){
          if(groupsize == 0)
            startGr = i;
          ids[groupsize] = J[i];
          groupsize++;
        }
      }
      if(groupsize > 1){
        //Sort indices within each group//
        base_arr = (double *) malloc(groupsize*sizeof(double));
        int indexSort[groupsize];
        for(i = 0; i < groupsize; i++){
          indexSort[i] = i;
          base_arr[i] = r[ids[i]];
        }
        qsort(indexSort, groupsize, sizeof(*indexSort), comparCTR);
        if(fabs(r[ids[indexSort[groupsize - 1]]]) > fabs(r[ids[indexSort[0]]]))
          reverseArrayCTR(indexSort, 0, groupsize-1);
        free(base_arr);
        //End of sort//

        R = 0; Dold = 0; Dnew = 0; Nold = 0; Nnew = 0; Nnewsq = 0; l = 0;
        memset(aNew, 0, p*sizeof(double)); memset(z, 0, n*sizeof(double));

        ////////// Split a given group //////////

        for (i = 0; i < groupsize; i++){
          t_j = ids[indexSort[i]];
          J[i + startGr] = t_j;
          index1 = n*t_j;
          qvalue = 0; vvalue = Dbestinv*v[t_j];
          for(j = 0; j < n; j++){
            index2 = index1 + j;
            zvalue = z[j];
            evalue = E[index2];
            evalue = evalue + vvalue*e_z[j];
            E[index2] = evalue;
            qvalue = qvalue + evalue*zvalue;
            zvalue = zvalue + X[index2];
            z[j] = zvalue;
          }
          Nnew = Nold + w[t_j];
          Dnew = Dold + 2*qvalue + u[t_j];
          Nnewsq = pow(Nnew, 2);
          if(Nnewsq > 1e-10)
            R = Nnewsq/Dnew;
          else
            R = 0;
          Nold = Nnew; Dold = Dnew; aNew[t_j] = 1; l = l + 1;
          if(R > Rbest){
            Rbest = R; Nbest = Nnew; Dbest = Dnew; splitGroup = k; splitPoint = l; splitsize = groupsize;
            memcpy(abest1, aNew, p * sizeof (double));
            memcpy(zbest, z, n * sizeof (double));
          }
        }
      }
      if(splitGroup == k){
        F77_CALL(daxpy)(dim, cdm1, abest1, ci1, aNew, ci1);
        memcpy(abest2, aNew, p * sizeof (double));
      }
    }

    //////////If SSE improves //////////
    if(Rbest > 0){
      memcpy(ZOld, Z, nngr*sizeof (double));
      zz = 0;
      index1 = iterNext * n;
      for(j = 0; j < n; j++){
        Z[index1 + j] = zbest[j];
        zz = zz + pow(zbest[j], 2);
      }
      //////////(ngr * 1 vector) Zz = t(ZOld) * z//////////
      F77_CALL(dgemv)("T", cn, cngr, cd1, ZOld, cn, zbest, ci1, cd0, Zz, ci1);
      //////////(ngr * 1 vector) ZZinvZz = t(ZZinv) * Zz//////////
      F77_CALL(dgemv)("T", cngr, cngr, cd1, ZZinv, cngr, Zz, ci1, cd0, ZZinvZz, ci1);
      ////////////Rank-1 Update//////////
      ZzZZinvZz = 0;
      for(j = 0; j < iterNext; j++)
        ZzZZinvZz = ZzZZinvZz + Zz[j]*ZZinvZz[j];
      denum = 1 / (zz - ZzZZinvZz);
      F77_CALL(dger)(cngr, cngr, &denum, ZZinvZz, ci1, ZZinvZz, ci1, ZZinv, cngr);
      index1 = iterNext*ngr;
      for (i = 0; i < iterNext; i++) {
        invvalue = -denum*ZZinvZz[i];
        ZZinv[iterNext + i*ngr] = invvalue;
        ZZinv[index1 + i] = invvalue;
      }
      ZZinv[index1 + iterNext] = denum;
      ////////////////////////////////////////
      memcpy(e_z, zbest, n*sizeof (double));
      Dbestinv = -1/Dbest;
      //////////e_z = zbest - ZtrOld * ZZinvZz//////////
      F77_CALL(dgemv)("N", cn, cngr, cdm1, ZOld, cn, ZZinvZz, ci1, cd1, e_z, ci1);
      //////////v = t(E) * e_z//////////
      F77_CALL(dgemv)("T", cn, dim,  cd1, E, cn, e_z, ci1, cd0, v, ci1);

      ////////// Update the matrices for search //////////
      abestXY = 0; l = 0;
      for (i = 0; i < p; i++){
        vvalue = v[i]; wvalue = w[i]; uvalue = u[i];
        wvalue = wvalue + Dbestinv*Nbest*vvalue;
        uvalue = uvalue + Dbestinv*pow(vvalue, 2);
        A[iter*p + i] = abest1[i];
        abestXY = abestXY + abest1[i]*XY[i];
        GroupTree[iter*p + i] = Groups[i];
        JTree[iter*p + i] = J[i];
        if (Groups[i] == splitGroup){
          isSplit[iter] = splitGroup;
          if (l >= splitPoint){
            if (Groups[i] == nongr[0])
              nongr[0] = iterNext;

            Groups[i] = iterNext;
            GroupTree[iter*p + i] = iterNext;
          }else{
            if(Groups[i] == nongr[0]){
              if(splitsize == splitPoint){
                nongr[0] = iterNext;
              }
            }
          }
          l = l + 1;
        }
        if(uvalue == 0){
          r[i] = 0;
        }else{
          r[i] = copysign(1.0, wvalue)*pow(wvalue, 2)/uvalue;
        }
        u[i] = uvalue; w[i] = wvalue;
      }
      ZY[iterNext] = abestXY;

      ////////// Compute SSE //////////
      F77_CALL(dgemv)("T", cngr, cngr, cd1, ZZinv, cngr, ZY, ci1, cd0, ZZinvZY, ci1);
      F77_CALL(dgemv)("N", cn, cngr, cd1, Z, cn, ZZinvZY, ci1, cd0, Yhat, ci1);
      F77_CALL(daxpy)(cn, cdm1, Y, ci1, Yhat, ci1);
      SSE = pow(F77_CALL(dnrm2)(cn, Yhat, ci1), 2);
      Rsq[iter] = 1 - SSE/SST;

      ////////// Update coefficients //////////
      memset(alpha, 0, ngroups*sizeof(double)); memset(beta, 0, p*sizeof(double));
      for (k = 0; k < iterNext; k++)
        alpha[k] = ZZinvZY[k + 1];
      F77_CALL(dgemv)("N", dim, cngroups, cd1, A, dim, alpha, ci1, cd0, beta, ci1);
      int pp1 = p + 1;
      betaTree[iter*pp1] = ZZinvZY[0];
      for(i = 1; i < pp1; i++)
        betaTree[iter*pp1 + i] = beta[i - 1];
    }
  }
  free(A); free(ZOld); free(ZZinv); free(Zz); free(ZZinvZz); free(ZY); free(XY); free(ZZinvZY); free(Yhat);
  free(aNew); free(abest1); free(abest2); free(zbest); free(z);
  free(v); free(E); free(w); free(u); free(r); free(Groups); free(J);
  free(colsum); free(alpha); free(beta);
}











