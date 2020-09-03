
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
static int compar (const void *a, const void *b){
  int aa = *((int *) a), bb = *((int *) b);
  if (base_arr[aa] < base_arr[bb])
    return -1;
  else if (base_arr[aa] == base_arr[bb])
    return 0;
  else
    return 1;
}
void reverseArray(int arr[], int start, int end){
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

void CTRwithCV(double *X, double *Y, int *dim, int *obs, int *foldid, int *nfolds, int *Kmax,
                   int *K, double *sumSqErr, int* Groupsf, int* Jf, double* Zf, double* ZZinvf){

  int i, j, k, f, l, iter, iterNext, index1, index2, index3, c1, c2;
  int ngroupsmax = *Kmax; int ngroupsmin = *K; int n = *obs; int p = *dim; int nfold = *nfolds; int one = 1; int ngr = ngroupsmax + 1;
  double d1 = 1.0; double d0 = 0.0; double dm1 = -1.0;
  const double *cd1 = &d1; const double *cd0 = &d0; const double *cdm1 = &dm1; const int *ci1 = &one; const int *cngr = &ngr;

  double *SSEcv = (double *) malloc(nfold*ngroupsmax*sizeof( double )); double *SSTcv = (double *) malloc(nfold*sizeof( double ));
  double *w = (double *) malloc(p*sizeof(double)); double *u = (double *) malloc(p*sizeof(double));
  double *r = (double *) malloc(p*sizeof(double)); double *v = (double *) malloc(p*sizeof(double));
  int *Groups = (int *) malloc(p*sizeof(int)); int *J = (int *) malloc(p*sizeof(int));
  double *ZZinv = (double *) malloc(ngr*ngr*sizeof( double )); double *Zz = (double *) malloc(ngr*sizeof( double ));
  double *ZZinvZz = (double *) malloc(ngr*sizeof( double )); double *ZZinvZY = (double *) malloc(ngr*sizeof( double ));
  double *XY = (double *) malloc(p*sizeof(double)); double *ZY = (double *) malloc(ngr*sizeof(double));
  double *abest = (double *) malloc(p*sizeof(double)); double *anew = (double *) malloc(p*sizeof(double));
  double *colsum = (double *) malloc(p*sizeof(double));

  double ZzZZinvZz, denum, invvalue, zzbest, abestXY;
  double Rbest, Nbest, Dbest, R, Dold, Dnew, Nold, Nnew, Nnewsq, Dbestinv = 0;
  double xvalue, yvalue, uvalue, wvalue, evalue, vvalue, qvalue, zvalue, xyvalue;
  int splitGroup, splitPoint, groupsize, t_j;
  double SSE, colmean, sumYtest, sumY, SSTtest;

  Rprintf("Start cross-validation\n");

  for (f = 0; f < nfold; f++) {
    int ntest = 0;
    for (j = 0; j < n; j++){
      if (foldid[j] == f)
        ntest = ntest + 1;
    }
    int ntr = n - ntest; const int *cntest = &ntest; const int *cntr = &ntr;
    int ntrp = p*ntr; int intrgr = ntr*ngr;

    double *E = (double *) malloc(ntrp*sizeof(double)); double *e_z = (double *) malloc(ntr*sizeof(double));
    double *Xtrain = (double *) malloc(ntrp*sizeof(double)); double *Xtest = (double *) malloc(ntest*p*sizeof(double));
    double *Ytest = (double *) malloc(ntest*sizeof(double)); double *Ytesthat = (double *) malloc(ntest*sizeof(double));
    double *Ytrain = (double *) malloc(ntr*sizeof(double)); double *ZtrOld = (double *) malloc(intrgr*sizeof(double));
    double *Ztr = (double *) malloc(intrgr*sizeof(double)); double *Ztest = (double *) malloc(ntest*ngr*sizeof(double));
    double *ztr = (double *) malloc(ntr*sizeof(double)); double *ztest = (double *) malloc(ntest*sizeof(double));
    double *zbest = (double *) malloc(ntr*sizeof(double)); double *zbesttest = (double *) malloc(ntest*sizeof(double));

    memset(E, 0, ntrp*sizeof(double)); memset(e_z, 0, ntr*sizeof(double));
    memset(Xtrain, 0, ntrp*sizeof(double)); memset(Xtest, 0, ntest*p*sizeof(double));
    memset(Ytest, 0, ntest*sizeof(double)); memset(Ytesthat, 0, ntest*sizeof(double));
    memset(Ytrain, 0, ntr*sizeof(double)); memset(ZtrOld, 0, intrgr*sizeof(double));
    memset(Ztr, 0, intrgr*sizeof(double)); memset(Ztest, 0, ntest*ngr*sizeof(double));
    memset(ZZinv, 0, ngr*ngr*sizeof(double));
    memset(Zz, 0, ngr*sizeof(double)); memset(ZZinvZz, 0, ngr*sizeof(double));
    memset(ZY, 0, ngr*sizeof(double)); memset(ZZinvZY, 0, ngr*sizeof(double));

    ZZinv[0] = (double) 1/ntr;
    sumYtest = 0; sumY = 0; SSTtest = 0;

    for (i = 0; i < p; i++){
      colsum[i] = 0; colmean = 0;
      c1 = 0; c2 = 0;
      index1 = i*n; index2 = i*ntr; index3 = i*ntest;
      for (j = 0; j < n; j++){
        if (foldid[j] == f){
          Xtest[c1 + index3] = X[j + index1];
          c1++;
        }else{
          xvalue = X[j + index1];
          Xtrain[c2 + index2] = xvalue;
          colmean = colmean + xvalue;
          c2++;
        }
      }
      colsum[i] = colmean/ntr;
    }

    c1 = 0; c2 = 0;
    for (j = 0; j < n; j++){
      if (foldid[j] == f){
        yvalue = Y[j];
        Ztest[c1] = 1; Ytest[c1] = yvalue; sumYtest = sumYtest + yvalue;
        Zf[f*n*ngr + j] = 1;
        c1++;
      }else{
        yvalue = Y[j];
        Ztr[c2] = 1; Ytrain[c2] = yvalue; sumY = sumY + yvalue;
        Zf[f*n*ngr + j] = 1;
        c2++;
      }
    }
    ZY[0] = sumY;
    sumYtest = (double) sumYtest/ntest;
    for (j = 0; j < ntest; j++)
      SSTtest = SSTtest + pow(Ytest[j] - sumYtest, 2);
    SSTcv[f] = SSTtest;

    if(ngroupsmin > 0){
      int index4 = f*n*ngr;
      for (i = 1; i < ngr; i++){
        c1 = 0; c2 = 0;
        index1 = i*n; index2 = i*ntr; index3 = i*ntest;
        for (j = 0; j < n; j++){
          if (foldid[j] == f){
            Ztest[index3 + c1] = Zf[index4 + index1 + j];
            c1++;
          }else{
            Ztr[index2 + c2] = Zf[index4 + index1 + j];
            c2++;
          }
        }
      }
      F77_CALL(dgemv)("T", cntr, cngr, cd1, Ztr, cntr, Ytrain, ci1, cd0, ZY, ci1);
      index4 = f*ngr*ngr;
      for (i = 0; i < ngr; i++)
        for (j = 0; j < ngr; j++)
          ZZinv[i*ngr + j] = ZZinvf[index4 + i*ngr + j];

      double *ZX = (double *) malloc(ngr*p*sizeof(double)); double *ZZinvZX = (double *) malloc(ngr*ntr*sizeof(double));
      memset(ZX, 0, ngr*p*sizeof(double)); memset(ZZinvZX, 0, ngr*ntr*sizeof(double));

      F77_CALL(dgemm)("T", "N", cngr, dim, cntr, cd1, Ztr, cntr, Xtrain, cntr, cd0, ZX, cngr);
      F77_CALL(dgemm)("T", "N", cngr, dim, cngr, cd1, ZZinv, cngr, ZX, cngr, cd0, ZZinvZX, cngr);
      F77_CALL(dgemm)("N", "N", cntr, dim, cngr, cd1, Ztr, cntr, ZZinvZX, cngr, cd0, E, cntr);

      free(ZX); free(ZZinvZX);
    }

    if(ngroupsmin == 0){
      for (i = 0; i < p; i++){
        uvalue = 0; wvalue = 0; xyvalue = 0;
        index1 = i*ntr;
        Groups[i] = 0; J[i] = i;
        u[i] = 0; w[i] = 0; XY[i] = 0;
        colmean = colsum[i];
        for(j = 0; j < ntr; j++){
          index2 = index1 + j;
          xvalue = Xtrain[index2];
          yvalue = Ytrain[j];
          evalue = xvalue - colmean;
          E[index2] = evalue;
          xyvalue = xyvalue + xvalue*yvalue;
          uvalue = uvalue + evalue*evalue;
          wvalue = wvalue + evalue*yvalue;
        }
        u[i] = uvalue; w[i] = wvalue; XY[i] = xyvalue;
        if(u[i] == 0){
          r[i] = 0;
        }else{
          r[i] = copysign(1.0, w[i])*pow(w[i], 2)/u[i];
        }
      }
    }else{
      for (i = 0; i < p; i++){
        uvalue = 0; wvalue = 0; xyvalue = 0;
        index1 = i*ntr;
        Groups[i] = Groupsf[f*p + i]; J[i] = Jf[f*p + i];
        u[i] = 0; w[i] = 0; XY[i] = 0;
        colmean = colsum[i];
        for(j = 0; j < ntr; j++){
          index2 = index1 + j;
          xvalue = Xtrain[index2];
          yvalue = Ytrain[j];
          E[index2] = xvalue - E[index2];
          evalue = E[index2];
          xyvalue = xyvalue + xvalue*yvalue;
          uvalue = uvalue + evalue*evalue;
          wvalue = wvalue + evalue*yvalue;
        }
        u[i] = uvalue; w[i] = wvalue; XY[i] = xyvalue;
        if(u[i] == 0){
          r[i] = 0;
        }else{
          r[i] = copysign(1.0, w[i])*pow(w[i], 2)/u[i];
        }
      }
    }

    for (iter = ngroupsmin; iter < ngroupsmax; iter++) {
      iterNext = iter + 1;
      Rbest = 0; Nbest = 0; Dbest = 0; splitGroup = 0; splitPoint = 0;
      memset(abest, 0, p*sizeof(double)); memset(zbest, 0, ntr*sizeof(double)); memset(zbesttest, 0, ntest*sizeof(double));

      //Find the split of each group//
      for (k = 0; k < iterNext; k++){

        //Determine the indices of the predictors in each group//
        int ids[p]; groupsize = 0; int startGr = 0;
        for(i = 0; i < p; i++){
          if(Groups[i] == k){
            if(groupsize == 0)
              startGr = i;
            ids[groupsize] = J[i];
            groupsize++;
          }
        }
        //Indices of the predictors are known//
        if(groupsize > 1){
          //Sort indices within each group//
          base_arr = (double *) malloc(groupsize*sizeof(double));
          int sortedindices[groupsize];
          for(i = 0; i < groupsize; i++){
            sortedindices[i] = i;
            base_arr[i] = r[ids[i]];
          }
          qsort(sortedindices, groupsize, sizeof(*sortedindices), compar);
          if(fabs(r[ids[sortedindices[groupsize - 1]]]) > fabs(r[ids[sortedindices[0]]]))
            reverseArray(sortedindices, 0, groupsize-1);
          free(base_arr);
          //End of sorted//

          memset(ztr, 0, ntr*sizeof(double)); memset(ztest, 0, ntest*sizeof(double)); memset(anew, 0, p*sizeof(double));
          R = 0; Dold = 0; Dnew = 0; Nold = 0; Nnew = 0; Nnewsq = 0; l = 0;

          for (i = 0; i < groupsize; i++){
            t_j = ids[sortedindices[i]];
            J[i + startGr] = t_j;
            qvalue = 0; vvalue = Dbestinv*v[t_j];
            index1 = ntr*t_j; index2 = ntest*t_j;
            for(j = 0; j < ntr; j++){
              index3 = j + index1;
              evalue = E[index3];
              evalue = evalue + vvalue*e_z[j];
              E[index3] = evalue;
              xvalue = Xtrain[index3];
              zvalue = ztr[j];
              qvalue = qvalue + zvalue*evalue;
              zvalue = zvalue + xvalue;
              ztr[j] = zvalue;
            }
            for(j = 0; j < ntest; j++){
              ztest[j] = ztest[j] + Xtest[j + index2];
            }
            Nnew = Nold + w[t_j];
            Dnew = Dold + 2*qvalue + u[t_j];
            Nnewsq = pow(Nnew, 2);
            if(Nnewsq > 1e-10)
              R = Nnewsq/Dnew;
            else
              R = 0;

            Nold = Nnew; Dold = Dnew; l = l + 1; anew[t_j] = 1;
            if(R > Rbest){
              Rbest = R; Nbest = Nnew; Dbest = Dnew; splitGroup = k; splitPoint = l;
              memcpy(abest, anew, p*sizeof (double));
              memcpy(zbest, ztr, ntr*sizeof(double));
              memcpy(zbesttest, ztest, ntest*sizeof(double));
            }
          }
        }
      }

      if(Rbest > 0){
        memcpy(ZtrOld, Ztr, intrgr*sizeof (double));
        index1 = ntr*iterNext;
        index2 = ntest*iterNext;
        zzbest = 0;
        c1 = 0; c2 = 0;
        for(j = 0; j < n; j++){
          if (foldid[j] == f){
            Ztest[index2 + c1] = zbesttest[c1];
            Zf[f*n*ngr + n*iterNext + j] = zbesttest[c1];
            c1++;
          }else{
            zvalue = zbest[c2];
            Ztr[index1 + c2] = zvalue;
            zzbest = zzbest + pow(zvalue, 2);
            Zf[f*n*ngr + n*iterNext + j] = zvalue;
            c2++;
          }
        }
        F77_CALL(dgemv)("T", cntr, cngr, cd1, ZtrOld, cntr, zbest, ci1, cd0, Zz, ci1);
        F77_CALL(dgemv)("T", cngr, cngr, cd1, ZZinv, cngr, Zz, ci1, cd0, ZZinvZz, ci1);
        //////////////////////////////////Calculate inverse//////////////////////////////////
        ZzZZinvZz = 0;
        for(j = 0; j < iterNext; j++)
          ZzZZinvZz = ZzZZinvZz + Zz[j]*ZZinvZz[j];

        denum = 1 / (zzbest - ZzZZinvZz);

        F77_CALL(dger)(cngr, cngr, &denum, ZZinvZz, ci1, ZZinvZz, ci1, ZZinv, cngr);
        index1 = iterNext*ngr;
        for (i = 0; i < iterNext; i++) {
          invvalue = -denum*ZZinvZz[i];
          ZZinv[iterNext + i*ngr] = invvalue;
          ZZinv[index1 + i] = invvalue;
        }
        ZZinv[index1 + iterNext] = denum;

        /////////////////////////////////////////////////////////////////////////////////////
        memcpy(e_z, zbest, ntr*sizeof (double));
        Dbestinv = -1/Dbest;
        // // //e_z = zbest - ZtrOld * ZZinvZz// // //
        F77_CALL(dgemv)("N", cntr, cngr, cdm1, ZtrOld, cntr, ZZinvZz, ci1, cd1, e_z, ci1);
        // // //v = t(E) * e_z// // //
        F77_CALL(dgemv)("T", cntr, dim, cd1, E, cntr, e_z, ci1, cd0, v, ci1);
        ///////////////Update Matrices///////////////////////////
        abestXY = 0; l = 0;

        for (i = 0; i < p; i++){
          vvalue = v[i]; uvalue = u[i]; wvalue = w[i];
          wvalue = wvalue + Dbestinv*Nbest*vvalue;
          uvalue = uvalue + Dbestinv*pow(vvalue, 2);
          abestXY = abestXY + abest[i]*XY[i];
          if(uvalue == 0){
            r[i] = 0;
          }else{
            r[i] = copysign(1.0, wvalue)*pow(wvalue, 2)/uvalue;
          }
          u[i] = uvalue;
          w[i] = wvalue;
          if(Groups[i] == splitGroup){
            if (l >= splitPoint)
              Groups[i] = iterNext;
            l = l + 1;
          }
        }
        ZY[iterNext] = abestXY;
        ////////////////////////////////////////////////////////////////////
        F77_CALL(dgemv)("T", cngr, cngr, cd1, ZZinv, cngr, ZY, ci1, cd0, ZZinvZY, ci1);
        F77_CALL(dgemv)("N", cntest, cngr, cd1, Ztest, cntest, ZZinvZY, ci1, cd0, Ytesthat, ci1);
        F77_CALL(daxpy)(cntest, cdm1, Ytest, ci1, Ytesthat, ci1);
        SSE = F77_CALL(dnrm2)(cntest, Ytesthat, ci1);
        ////////////////////////////////////////////////////////////////////
        SSEcv[f + nfold*iter] = pow(SSE, 2);

      }else{
        SSEcv[f + nfold*iter] = SSEcv[f + nfold*(iter-1)];
      }
    }
    free(Xtrain); free(Xtest); free(Ytesthat); free(Ytest); free(Ytrain); free(E); free(Ztr); free(Ztest); free(ZtrOld);
    free(ztr); free(ztest); free(zbest); free(zbesttest); free(e_z);
    for (i = 0; i < p; i++){
      Groupsf[f*p + i] = Groups[i];
      Jf[f*p + i] = J[i];
    }
    for (i = 0; i < ngr; i++){
      for (j = 0; j < ngr; j++){
        ZZinvf[f*ngr*ngr + i*ngr + j] = ZZinv[i*ngr + j];
      }
    }

  }
  free(w); free(u); free(r); free(v); free(Groups); free(J); free(ZZinv); free(Zz); free(ZZinvZz); free(ZZinvZY); free(XY); free(ZY);
  free(abest); free(anew); free(colsum);

  Rprintf("End cross-validation\n");
  Rprintf("...........................\n");

  double sumSSEcv = 0;
  for (i = ngroupsmin; i < ngroupsmax; i++){
    sumSSEcv = 0;
    for(j = 0; j < nfold; j++)
      sumSSEcv = sumSSEcv + SSEcv[j + i*nfold]/SSTcv[j];
    sumSqErr[i] = sumSSEcv/nfold;
  }

  double minCV = 1.0 /0.0;
  int bestSize = 0;
  for (k = 0; k < ngroupsmax; k++){

    if(sumSqErr[k] < minCV){
      bestSize = k + 1;
      minCV = sumSqErr[k];
    }
  }
  K[0] = bestSize;
  Rprintf("Optimal group number : %i\n", K[0]);
  Rprintf("...........................\n");
  free(SSEcv); free(SSTcv);
}

