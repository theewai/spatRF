#include "R.h"
#include "Rinternals.h"
#include "R_ext/Lapack.h"

//void F77_NAME(llc_f)(double *x, int n, double *l, double *a, double *ret);

// Finds x^T A^{- 1} x using DPOTRS to avoid inverse
extern "C"
void matsmash(double *a, int *n, double *x, double *result,int *index){
    int info;
    int ione = 1;
    double *savex = (double *) R_alloc(n[0], sizeof(double));
    F77_CALL(dpotrf)("L", n, a, n, &info);
    if (info != 0) {
        result[*index]=-1;
        return;
        //error("Cholesky decomposition failed");
    }
    F77_CALL(dcopy)(n, x, &ione, savex, &ione);
    F77_CALL(dpotrs)("L", n, &ione, a, n, x, n, &info);
    if (info != 0) {
        result[*index]=-1;
        return;
        //error("Cholesky decomposition failed");
    }
    result[*index] = F77_CALL(ddot)(n, x, &ione, savex, &ione);
}

extern "C"
SEXP C_splits(SEXP x, SEXP numperleaf,SEXP splitleaf , SEXP leafmem, SEXP smpl, SEXP r,SEXP sortind, SEXP kern,SEXP omegaY, SEXP L,SEXP numpergroup, SEXP group){
int k = length(splitleaf);
int n = length(leafmem);
int s = length(smpl);
int nummaxsplit = s* n;
int numsplits=0;
int strtind[k];
int obsind,covarnum,obsleaf,xnum,leafind,xnumnext,yind,yindold;
int leaforder,groupind,groupindold;
int count[k];
//int strtj = 0;
//int lmadjust = 1;
int num_groups = length(numpergroup)/k;
int copyn,copyns;
int max_small, min_large,ci;
int count_group[num_groups], left_group[num_groups];
double pen[num_groups*num_groups],copypen[num_groups*num_groups];
double weights[num_groups],copyweights[num_groups];

int *pnpl, *psl,*plm,*psmpl,*pr,*pcovar,*pleaf,*psortind,*pcli,*pns,*pgroup,*pnpg,*pdebug,*pdebugsec;
double *pctpt,*px,*ploss, *pkern,*poy,*pL;

pnpl = INTEGER(numperleaf);
pnpg = INTEGER(numpergroup);
psl = INTEGER(splitleaf);
plm = INTEGER(leafmem);
psmpl = INTEGER(smpl);
pr = INTEGER(r);
pgroup = INTEGER(group);
px = REAL(x);
psortind = INTEGER(sortind);
pkern = REAL(kern);
poy = REAL(omegaY);
pL = REAL(L);

SEXP ctpt = PROTECT(allocVector(REALSXP, nummaxsplit));
SEXP covar = PROTECT(allocVector(INTSXP, nummaxsplit));
SEXP leaf = PROTECT(allocVector(INTSXP, nummaxsplit));
SEXP covleafind = PROTECT(allocVector(INTSXP, n));
SEXP vec = PROTECT(allocVector(VECSXP, 7));
SEXP ns = PROTECT(allocVector(INTSXP, 1));
SEXP loss = PROTECT(allocVector(REALSXP, nummaxsplit));
SEXP debug =PROTECT(allocVector(INTSXP,nummaxsplit));
SEXP debugsec =PROTECT(allocVector(INTSXP,nummaxsplit));

pcli = INTEGER(covleafind);
pctpt = REAL(ctpt);
pcovar = INTEGER(covar);
pleaf = INTEGER(leaf);
pns = INTEGER(ns);
ploss = REAL(loss);
pdebug = INTEGER(debug);
pdebugsec = INTEGER(debugsec);

// Set up counter indices for each leaf
strtind[0] = 0;
for(int t = 1; t < k; t++){
  strtind[t] = strtind[t-1] + pnpl[t-1];
}

// For Each Covariate in the sample
for (int i = 0; i < s; i++) {
  // off by one since R starts at 1
  covarnum = psmpl[i]-1;

  for(int d=0; d<k; d++){
    count[d] = 0;
  }

  // Sort covariate first by leaf, then in order of value
  for(int a=0; a<n; a++){
    // sortind off by one since R starts at 1
    obsind = psortind[a + n * covarnum] - 1;
    if(plm[obsind]!=0){
        obsleaf = plm[obsind] - 1;
        leaforder = strtind[obsleaf] + count[obsleaf];
        pcli[leaforder] = obsind;
        count[obsleaf] += 1;
    }
  }
  // covleafind now contains the order of observations,
  // first sorted by leaf, then by covariate

// Initialize the leaf counter
  for(int j = 0; j < k; j++){
    
    ci = 0;
    
    // Set the size of the smallest left leaf to zero
    max_small=0;
    // Find the size of the largest group
    min_large=n;

    // Reset Weights and Penalty matrix to the original form
    for(int ri = 0; ri < num_groups; ri++){
      weights[ri]=0;
      count_group[ri]=0;
      left_group[ri] = pnpg[num_groups*(psl[j]-1)+ri];
      pen[ri * num_groups + ri] = pL[ri * num_groups + ri];
      if(left_group[ri] < min_large){
        min_large = left_group[ri];
      }
      for(int rii =0; rii < ri; rii++){
          pen[ri * num_groups + rii] = pen[rii * num_groups + ri] = pL[ri*num_groups + rii];
      }
    }

    // For each observation in the leaf
    while( min_large >= *pr ){

      // Get index of which observation is next lowest
      leafind = strtind[j] + ci;
      yind = pcli[leafind];
      groupind = pgroup[yind] - 1;
      count_group[groupind] += 1;
      left_group[groupind] -= 1;

      if(count_group[groupind] > max_small) {
        max_small = count_group[groupind];
      }

      if(left_group[groupind] < min_large) {
        min_large = left_group[groupind];
      }

      // Update OmegaY and (C^T Omega C +L)

      weights[groupind] += poy[yind];

      if(ci > 0){
          for(int bi = 0; bi < ci; bi++){
              yindold = pcli[ strtind[j] + bi ];
              groupindold = pgroup[yindold]-1;
              pen[groupind * num_groups + groupindold]= pen[groupindold * num_groups + groupind] +=  pkern[n * yindold + yind];
          }
      }

      pen[groupind * num_groups + groupind] += pkern[n * yind + yind];

      if(max_small >= *pr & min_large >= *pr){

        xnum = n*covarnum + pcli[leafind];
        xnumnext = n*covarnum + pcli[leafind +1];

        if(px[xnum] - px[xnumnext] < 0){
          pctpt[numsplits] = (px[xnum] + px[xnumnext])/2;
          pleaf[numsplits] = psl[j];
          pcovar[numsplits] = covarnum + 1;
          pdebug[numsplits]= max_small;
          pdebugsec[numsplits] = min_large;

          for(int copyi=0; copyi < num_groups;copyi++){
            copyweights[copyi] = weights[copyi];
            copypen[copyi * num_groups + copyi] = pen[copyi * num_groups + copyi];

            for(int copyj=0; copyj < copyi; copyj++){
              copypen[copyi * num_groups + copyj] = copypen[copyj * num_groups + copyi] = pen[copyi*num_groups + copyj];
            }
          }

          copyn = num_groups;
          copyns=numsplits;

          matsmash(copypen, &copyn, copyweights, ploss,&copyns);

          numsplits += 1;

        }

      }
        
      ci++;

      // End loop for each observation in the leaf
      }

    // End loop for each leaf
    }

  // End loop for each covariate
  }

pns[0] = numsplits;

SET_VECTOR_ELT(vec, 0, ctpt);
SET_VECTOR_ELT(vec, 1, covar);
SET_VECTOR_ELT(vec, 2, leaf);
SET_VECTOR_ELT(vec, 3, ns);
SET_VECTOR_ELT(vec, 4, loss);
SET_VECTOR_ELT(vec, 5, debug);
SET_VECTOR_ELT(vec, 6, debugsec);
UNPROTECT(9);

return vec;

}

extern "C"
SEXP C_splits1( SEXP x, SEXP numperleaf, SEXP splitleaf, SEXP leafmem,
               SEXP smpl,SEXP r, SEXP sortind, SEXP kern, SEXP omegaY){
  int k = length(splitleaf);
  int n = length(leafmem);
  int s = length(smpl);
  int numsplits=0;
  int strtind[k];
  int obsind,covarnum,leafnum,obsleaf,xnum,leafind,xnumnext,ynum,yind,yindold;
  int nummaxsplit,yold,leaforder;
  int count[k];
  int strtj = 0;
  int lmadjust = 1;
  double runningnum, runningden;
  
  int *pnpl, *psl,*plm,*psmpl,*pr,*pcovar,*pleaf,*psortind,*pcli,*pns;
  double *pctpt,*px,*plossnum,*plossden,*pkern,*poy ;
  
  pnpl = INTEGER(numperleaf);
  psl = INTEGER(splitleaf);
  plm = INTEGER(leafmem);
  psmpl = INTEGER(smpl);
  pr = INTEGER(r);
  px = REAL(x);
  psortind = INTEGER(sortind);
  pkern = REAL(kern);
  poy = REAL(omegaY);
  
  nummaxsplit = s* n; 
  SEXP ctpt = PROTECT(allocVector(REALSXP, nummaxsplit));
  SEXP covar = PROTECT(allocVector(INTSXP, nummaxsplit));
  SEXP leaf = PROTECT(allocVector(INTSXP, nummaxsplit));
  SEXP covleafind = PROTECT(allocVector(INTSXP, n));
  SEXP vec = PROTECT(allocVector(VECSXP, 5));
  SEXP ns = PROTECT(allocVector(INTSXP, 1));
  SEXP lossnum = PROTECT(allocVector(REALSXP, nummaxsplit));
  SEXP lossden = PROTECT(allocVector(REALSXP, nummaxsplit));
  
  pcli = INTEGER(covleafind);
  pctpt = REAL(ctpt);
  pcovar = INTEGER(covar);
  pleaf = INTEGER(leaf);
  pns = INTEGER(ns);
  plossnum = REAL(lossnum);
  plossden = REAL(lossden);
  
  // Set up counter indices for each leaf
  strtind[0] = 0;
  for(int t = 1; t < k; t++){
    strtind[t] = strtind[t-1] + pnpl[t-1];
  }
  
  // Check if there are leaves that shouldn't be split
  if(psl[0] ==0){
    strtj =1;
    lmadjust = 0;
  } 
  
  // For Each Covariate in the sample
  for ( int i = 0; i < s; i++ ) {
    covarnum = psmpl[i]-1;    
    for(int d=0; d<k; d++){
      count[d] = 0;
    }
    for(int check =0; check < n; check++){
      pcli[check]=check;
    }
    
    // Sort covariate first by leaf, then in order of value
    for( int a=0; a < n; a++){
      obsind = psortind[ a + n * covarnum ]-1;
      obsleaf = plm[ obsind] - lmadjust;
      leaforder = strtind[ obsleaf ] + count[ obsleaf ];
      pcli[leaforder] = obsind;
      count[obsleaf]+=1;
    }
    // covleafind now contains the order of observations, first sorted by leaf, then by covariate
    
    // Initialize the leaf counter
    for( int j = strtj; j < k; j++ ){
      leafnum = psl[j];
      runningnum = 0;
      runningden = 0;
      
      if( (*pr-1) > 0 ){
        for(int ci = 0; ci < (*pr-1); ci++){
          yind = pcli[ strtind[j] + ci ];
          runningnum += poy[ yind ];
          if( ci > 0){
            for( int cii = 0; cii < ci; cii++){
              yindold = pcli[ strtind[j] + cii ];
              runningden += 2 * pkern[n * yind + yindold];
            }
          }
          runningden += pkern[n * yind + yind];
          
        }
      }
      
      // For each observation in the leaf
      for( int b = *pr-1; b < (pnpl[j] - *pr); b++){
        
        // Get index of which observation is next lowest               
        leafind = strtind[j] + b;
        ynum = pcli[ leafind ];
        xnum = n*covarnum + pcli[ leafind ];
        xnumnext = n*covarnum + pcli[ leafind +1 ];
        
        // Update running numerator and denominator total  
        runningnum += poy[ynum];
        
        if(b > 0){
          for( int bi = 0; bi < b; bi++){
            yold = pcli[ strtind[j] + bi ];
            runningden += 2 * pkern[n * ynum + yold];
          }
        }
        runningden +=  pkern[n * ynum + ynum];
        
        if(px[ xnum ] - px[ xnumnext] < 0){
          pctpt[numsplits] = (px[ xnum ] + px[ xnumnext])/2;
          pleaf[numsplits] = leafnum;
          pcovar[numsplits] = covarnum+1;
          plossnum[numsplits] = runningnum*runningnum/runningden;
          plossden[numsplits] = runningden;
          numsplits += 1;
        }
        
        // End loop for each observation in the leaf
      } 
      
      // End loop for each leaf  
    }
    
    // End loop for each covariate
  } 
  
  pns[0] = numsplits;
  SET_VECTOR_ELT(vec, 0, ctpt);
  SET_VECTOR_ELT(vec, 1, covar);
  SET_VECTOR_ELT(vec, 2, leaf);
  SET_VECTOR_ELT(vec, 3, ns);    
  SET_VECTOR_ELT(vec, 4, lossnum);
  UNPROTECT(8);
  
  return vec;
}
