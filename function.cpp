// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;
using namespace Eigen;

// identify the parental origin of off-spring alleles based on haplotype
// [[Rcpp::export]]
List distinguishC(MatrixXd gmm, MatrixXd gcc, int loci, MatrixXd hap) {
  int l = gmm.rows();
  int K = hap.rows();
  const VectorXd gm = gmm.col(loci - 1);
  const VectorXd gc = gcc.col(loci - 1);
  VectorXd gcm(l);
  gcm.fill(3);
  VectorXd gcp(l);
  gcp.fill(0);
  VectorXd phi(l);
  phi.fill(0);
  
  for(int u = 0; u < l; u++) {
    if (gc[u]==0) {
      gcm[u] = 0;
      gcp[u] = 0;
      phi[u] = 1;
    } 
    else if (gm[u] == 0 && gc[u] == 1) {
      gcm[u] = 0;
      gcp[u] = 1;
      phi[u] = 1;
    }
    else if (gm[u] == 1 && gc[u] == 1) {
      for(int k = 0; k < K; k++) {
        VectorXd hma = hap.row(k);
        VectorXd gmmu = gmm.row(u);
        VectorXd hmb = gmmu - hma;
        int r = 0;
        for(int kk = 0; kk < K; kk++) {
          VectorXd hapkk = hap.row(kk);
          int r0 = (hmb == hapkk) ? 1 : 0;
          r += r0;
        }
        if (r == 1) {
          VectorXd hcm = hma;
          VectorXd gccu = gcc.row(u);
          VectorXd hcp = gccu - hcm;
          int t = 0;
          for(int kk = 0; kk < K; kk++) {
            VectorXd hapkk = hap.row(kk);
            int t0 = (hcp == hapkk) ? 1 : 0;
            t += t0;
          }
          if (t == 1 && gcm[u] != hcm[loci-1]) {
            gcm[u] = hcm[loci-1];
            gcp[u] = hcp[loci-1];
            phi[u] = (phi[u] == 0) ? 1 : (phi[u] + 1);
          }
        }
      }
    }
    else if (gm[u] == 2 && gc[u] == 1) {
      gcm[u] = 1;
      gcp[u] = 0;
      phi[u] = 1;
    }
    else {
      gcm[u] = 1;
      gcp[u] = 1;
      phi[u] = 1;
    }
  }
  
  return List::create(_["gm"] = gm, _["gc"] = gc, _["gcm"] = gcm, _["gcp"] = gcp, _["phi"] = phi);
}

// calculate minor allel frequence based on haplotype
int whichrow_newC(VectorXd hcm, MatrixXd hap) {
  int K = hap.rows();
  int m;
  for(int kk = 0; kk < K; kk++) {
    VectorXd hapkk = hap.row(kk);
    int r0 = (hcm == hapkk) ? 1 : 0;
    if (r0 == 1) {
      m = kk;
    }
  }
  return m;
}
// [[Rcpp::export]]
List distinguish0C(MatrixXd gmm, MatrixXd gcc, int loci, MatrixXd hap, VectorXd ppi) {
  int l = gmm.rows();
  int K = hap.rows();
  VectorXd gcm(l);
  gcm.fill(3);
  VectorXd Pm(l);
  Pm.fill(0);
  VectorXd Pp(l);
  Pp.fill(0);
  for(int u = 0; u < l; u++) {
    for(int k = 0; k < K; k++) {
      VectorXd hma = hap.row(k);
      VectorXd gmmu = gmm.row(u);
      VectorXd hmb = gmmu - hma;
      int r = 0;
      for(int kk = 0; kk < K; kk++) {
        VectorXd hapkk = hap.row(kk);
        int r0 = (hmb == hapkk) ? 1 : 0;
        r += r0;
      }
      if (r == 1) {
        VectorXd hcm = hma;
        VectorXd gccu = gcc.row(u);
        VectorXd hcp = gccu - hcm;
        int t = 0;
        for(int kk = 0; kk < K; kk++) {
          VectorXd hapkk = hap.row(kk);
          int t0 = (hcp == hapkk) ? 1 : 0;
          t += t0;
        }
        if (t == 1) {
          gcm[u] = hcm[loci-1];
          if (gcm[u] == 1) {
            Pm[u] = Pm[u] + ppi[whichrow_newC(hma, hap)] * ppi[whichrow_newC(hmb, hap)] * ppi[whichrow_newC(hcp, hap)];
          }
          else {
            Pp[u] = Pp[u] + ppi[whichrow_newC(hma, hap)] * ppi[whichrow_newC(hmb, hap)] * ppi[whichrow_newC(hcp, hap)];
          }
        }
      }
    }
  }
  return List::create(_["Pm"] = Pm, _["Pp"] = Pp);
}

// identify the parental origin of off-spring alleles based on single SNP
// [[Rcpp::export]]
List distinguish_snpC(VectorXd gm, VectorXd gc) {
  int n = gm.size();
  VectorXd gcm(n);
  gcm.fill(3);
  VectorXd gcp(n);
  gcp.fill(0);
  VectorXd phi(n);
  phi.fill(0);
  
  for(int u = 0; u < n; u++) {
    if (gc[u]==0) {
      gcm[u] = 0;
      gcp[u] = 0;
      phi[u] = 1;
    } 
    else if (gm[u] == 0 && gc[u] == 1) {
      gcm[u] = 0;
      gcp[u] = 1;
      phi[u] = 1;
    }
    else if (gm[u] == 1 && gc[u] == 1) {
      gcm[u] = 3;
      gcp[u] = 3;
      phi[u] = 2;
    }
    else if (gm[u] == 2 && gc[u] == 1) {
      gcm[u] = 1;
      gcp[u] = 0;
      phi[u] = 1;
    }
    else {
      gcm[u] = 1;
      gcp[u] = 1;
      phi[u] = 1;
    }
  }
  return List::create(_["gm"] = gm, _["gc"] = gc, _["gcm"] = gcm, _["gcp"] = gcp, _["phi"] = phi);
}

// Pr(gc|gm;pa)
MatrixXd gc_on_gmC(double pa) {
  Matrix3d g;
  g(0,0) = 1. - pa;
  g(0,1) = pa;
  g(0,2) = 0;
  g(1,0) = (1. - pa)/2;
  g(1,1) = 1./2;
  g(1,2) = pa/2;
  g(2,0) = 0;
  g(2,1) = 1. - pa;
  g(2,2) = pa;
  return g;
}
// d(Pr(gc|gm;pa))|d(pa)
MatrixXd gc_on_gm_grC(double pa) {
  Matrix3d g;
  g(0,0) = -1.;
  g(0,1) = 1.;
  g(0,2) = 0;
  g(1,0) = -1./2;
  g(1,1) = 0;
  g(1,2) = 1./2;
  g(2,0) = 0;
  g(2,1) = -1.;
  g(2,2) = 1.;
  return g;
}

// design matrix test1
// [[Rcpp::export]]
MatrixXd design_mat1C(VectorXd gm, VectorXd gc, VectorXd gcm, VectorXd gcp, VectorXd x) {
  MatrixXd res(gm.size(),5);
  res.fill(1);
  for(int u = 0; u < gm.size(); u++) {
    res(u,1) = gm[u];
    res(u,2) = gc[u];
    res(u,3) = gcm[u] - gcp[u];
    res(u,4) = x[u];
  }
  return res;
}

MatrixXd design_mat2_1C(VectorXd gm, VectorXd gc, VectorXd gcm, VectorXd gcp, VectorXd x, VectorXd phi) {
  MatrixXd res(gm.size(),5);
  res.fill(1);
  for(int u = 0; u < gm.size(); u++) {
    if (phi[u] != 1) {
      gcm[u] = 1;
      gcp[u] = 0;
    }
  }
  for(int u = 0; u < gm.size(); u++) {
    res(u,1) = gm[u];
    res(u,2) = gc[u];
    res(u,3) = gcm[u] - gcp[u];
    res(u,4) = x[u];
  }
  return res;
}

MatrixXd design_mat2_2C(VectorXd gm, VectorXd gc, VectorXd gcm, VectorXd gcp, VectorXd x, VectorXd phi) {
  MatrixXd res(gm.size(),5);
  res.fill(1);
  for(int u = 0; u < gm.size(); u++) {
    if (phi[u] != 1) {
      gcm[u] = 0;
      gcp[u] = 1;
    }
  }
  for(int u = 0; u < gm.size(); u++) {
    res(u,1) = gm[u];
    res(u,2) = gc[u];
    res(u,3) = gcm[u] - gcp[u];
    res(u,4) = x[u];
  }
  return res;
}

// using haplotype to infer parent-of-origin
// Pr(Y=y|gmm,gcc,X;beta)
// [[Rcpp::export]]
VectorXd prYC(VectorXd y, MatrixXd gmm, MatrixXd gcc, VectorXd x, VectorXd beta, int m, MatrixXd hap, VectorXd ppi) {
  List F = distinguishC(gmm, gcc, m, hap);
  VectorXd gm = F["gm"];
  VectorXd gc = F["gc"];
  VectorXd gcm = F["gcm"];
  VectorXd gcp = F["gcp"];
  VectorXd phi = F["phi"];
  int n = y.size();
  int M = hap.cols();
  
  VectorXd Pm(n);
  Pm.fill(1);
  VectorXd Pp(n);
  Pp.fill(0);
  
  int r = 0;
  for (int u = 0; u < n; u++) {
    if (phi[u] != 1) {
      r +=1;
    }
  }
  MatrixXd gmm0(r,M);
  MatrixXd gcc0(r,M);
  
  MatrixXd Z = design_mat1C(gm, gc, gcm, gcp, x);
  MatrixXd Z1;
  MatrixXd Z2;
  
  if (phi.mean() == 1) {
    Z1 = Z;
    Z2 = Z;
  }
  else {
    int k = 0;
    for (int u = 0; u < n; u++) {
      if (phi[u] != 1) {
        for (int j = 0; j < M; j++) {
          gmm0(k,j) = gmm(u,j);
          gcc0(k,j) = gcc(u,j);
        }
        k += 1;
      }
    }
    
    Z1 = design_mat2_1C(gm, gc, gcm, gcp, x, phi);
    Z2 = design_mat2_2C(gm, gc, gcm, gcp, x, phi);
    
    List W = distinguish0C(gmm0, gcc0, m, hap, ppi);
    VectorXd Pm0 = W["Pm"];
    VectorXd Pp0 = W["Pp"];
    
    int t = 0;
    for (int u = 0; u < n; u++) {
      if (phi[u] != 1) {
        Pm[u] = Pm0[t];
        Pp[u] = Pp0[t];
        t += 1;
      }
    }
  }
  
  VectorXd P1(n);
  VectorXd P2(n);
  for (int u = 0; u < n; u++) {
    P1[u] = Pm[u] / (Pm[u] + Pp[u]);
    P2[u] = Pp[u] / (Pm[u] + Pp[u]);
  }
  
  VectorXd y1 = Z1 * beta;
  VectorXd y2 = Z2 * beta;
  
  VectorXd res(n);
  for (int u = 0; u < n; u++) {
    res[u] = P1[u] * exp(y[u]*y1[u]) / (1 + exp(y1[u])) + P2[u] * exp(y[u]*y2[u]) / (1 + exp(y2[u]));
  }
  return res;
}

MatrixXd prY_grC(VectorXd y, MatrixXd gmm, MatrixXd gcc, VectorXd x, VectorXd beta, int m, MatrixXd hap, VectorXd ppi) {
  List F = distinguishC(gmm, gcc, m, hap);
  VectorXd gm = F["gm"];
  VectorXd gc = F["gc"];
  VectorXd gcm = F["gcm"];
  VectorXd gcp = F["gcp"];
  VectorXd phi = F["phi"];
  int n = y.size();
  int M = hap.cols();
  
  VectorXd Pm(n);
  Pm.fill(1);
  VectorXd Pp(n);
  Pp.fill(0);
  
  int r = 0;
  for (int u = 0; u < n; u++) {
    if (phi[u] != 1) {
      r +=1;
    }
  }
  MatrixXd gmm0(r,M);
  MatrixXd gcc0(r,M);
  
  MatrixXd Z = design_mat1C(gm, gc, gcm, gcp, x);
  MatrixXd Z1;
  MatrixXd Z2;
  
  if (phi.mean() == 1) {
    Z1 = Z;
    Z2 = Z;
  }
  else {
    int k = 0;
    for (int u = 0; u < n; u++) {
      if (phi[u] != 1) {
        for (int j = 0; j < M; j++) {
          gmm0(k,j) = gmm(u,j);
          gcc0(k,j) = gcc(u,j);
        }
        k += 1;
      }
    }
    
    Z1 = design_mat2_1C(gm, gc, gcm, gcp, x, phi);
    Z2 = design_mat2_2C(gm, gc, gcm, gcp, x, phi);
    
    List W = distinguish0C(gmm0, gcc0, m, hap, ppi);
    VectorXd Pm0 = W["Pm"];
    VectorXd Pp0 = W["Pp"];
    
    int t = 0;
    for (int u = 0; u < n; u++) {
      if (phi[u] != 1) {
        Pm[u] = Pm0[t];
        Pp[u] = Pp0[t];
        t += 1;
      }
    }
  }
  
  VectorXd P1(n);
  VectorXd P2(n);
  for (int u = 0; u < n; u++) {
    P1[u] = Pm[u] / (Pm[u] + Pp[u]);
    P2[u] = Pp[u] / (Pm[u] + Pp[u]);
  }
  
  VectorXd y1 = Z1 * beta;
  VectorXd y2 = Z2 * beta;
  
  VectorXd betagr1(n);
  VectorXd betagr2(n);
  for (int u = 0; u < n; u++) {
    betagr1[u] = P1[u] * (2*y[u] - 1) * exp(y1[u]) / pow(1 + exp(y1[u]), 2);
    betagr2[u] = P2[u] * (2*y[u] - 1) * exp(y2[u]) / pow(1 + exp(y2[u]), 2);
  }
  
  MatrixXd betagr(n,5);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < 5; v++) {
      betagr(u,v) = Z1(u,v) * betagr1[u] + Z2(u,v) * betagr2[u];
    }
  }
  return betagr;
}

// [[Rcpp::export]]
MatrixXd log_prY_grC(VectorXd y, MatrixXd gmm, MatrixXd gcc, VectorXd x, VectorXd beta, int m, MatrixXd hap, VectorXd ppi) {
  List F = distinguishC(gmm, gcc, m, hap);
  VectorXd gm = F["gm"];
  VectorXd gc = F["gc"];
  VectorXd gcm = F["gcm"];
  VectorXd gcp = F["gcp"];
  VectorXd phi = F["phi"];
  int n = y.size();
  int M = hap.cols();
  
  VectorXd Pm(n);
  Pm.fill(1);
  VectorXd Pp(n);
  Pp.fill(0);
  
  int r = 0;
  for (int u = 0; u < n; u++) {
    if (phi[u] != 1) {
      r +=1;
    }
  }
  MatrixXd gmm0(r,M);
  MatrixXd gcc0(r,M);
  
  MatrixXd Z = design_mat1C(gm, gc, gcm, gcp, x);
  MatrixXd Z1;
  MatrixXd Z2;
  
  if (phi.mean() == 1) {
    Z1 = Z;
    Z2 = Z;
  }
  else {
    int k = 0;
    for (int u = 0; u < n; u++) {
      if (phi[u] != 1) {
        for (int j = 0; j < M; j++) {
          gmm0(k,j) = gmm(u,j);
          gcc0(k,j) = gcc(u,j);
        }
        k += 1;
      }
    }
    
    Z1 = design_mat2_1C(gm, gc, gcm, gcp, x, phi);
    Z2 = design_mat2_2C(gm, gc, gcm, gcp, x, phi);
    
    List W = distinguish0C(gmm0, gcc0, m, hap, ppi);
    VectorXd Pm0 = W["Pm"];
    VectorXd Pp0 = W["Pp"];
    
    int t = 0;
    for (int u = 0; u < n; u++) {
      if (phi[u] != 1) {
        Pm[u] = Pm0[t];
        Pp[u] = Pp0[t];
        t += 1;
      }
    }
  }
  
  VectorXd P1(n);
  VectorXd P2(n);
  for (int u = 0; u < n; u++) {
    P1[u] = Pm[u] / (Pm[u] + Pp[u]);
    P2[u] = Pp[u] / (Pm[u] + Pp[u]);
  }
  
  VectorXd y1 = Z1 * beta;
  VectorXd y2 = Z2 * beta;
  
  VectorXd res(n);
  for (int u = 0; u < n; u++) {
    res[u] = P1[u] * exp(y[u]*y1[u]) / (1 + exp(y1[u])) + P2[u] * exp(y[u]*y2[u]) / (1 + exp(y2[u]));
  }
  
  VectorXd betagr1(n);
  VectorXd betagr2(n);
  for (int u = 0; u < n; u++) {
    betagr1[u] = P1[u] * (2*y[u] - 1) * exp(y1[u]) / (pow(1 + exp(y1[u]), 2) * res[u]);
    betagr2[u] = P2[u] * (2*y[u] - 1) * exp(y2[u]) / (pow(1 + exp(y2[u]), 2) * res[u]);
  }
  
  MatrixXd betagr(n,5);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < 5; v++) {
      betagr(u,v) = Z1(u,v) * betagr1[u] + Z2(u,v) * betagr2[u];
    }
  }
  return betagr;
}

// using total probability function when both mother and child are heterozygous
// Pr(Y=y|gm,gc,X;beta,pa)
// [[Rcpp::export]]
VectorXd prY0C(VectorXd y, VectorXd gm, VectorXd gc, VectorXd x, VectorXd beta, double pa) {
  List F = distinguish_snpC(gm, gc);
  VectorXd gcm = F["gcm"];
  VectorXd gcp = F["gcp"];
  VectorXd phi = F["phi"];
  int n = y.size();
  
  MatrixXd Z = design_mat1C(gm, gc, gcm, gcp, x);
  MatrixXd Z1;
  MatrixXd Z2;
  
  if (phi.mean() == 1) {
    Z1 = Z;
    Z2 = Z;
  }
  else {
    Z1 = design_mat2_1C(gm, gc, gcm, gcp, x, phi);
    Z2 = design_mat2_2C(gm, gc, gcm, gcp, x, phi);
  }
  
  double P1 = 1 - pa;
  double P2 = pa;
  
  VectorXd y1 = Z1 * beta;
  VectorXd y2 = Z2 * beta;
  
  VectorXd res(n);
  for (int u = 0; u < n; u++) {
    res[u] = P1 * exp(y[u]*y1[u]) / (1 + exp(y1[u])) + P2 * exp(y[u]*y2[u]) / (1 + exp(y2[u]));
  }
  return res;
}

MatrixXd prY0_grC(VectorXd y, VectorXd gm, VectorXd gc, VectorXd x, VectorXd beta, double pa) {
  List F = distinguish_snpC(gm, gc);
  VectorXd gcm = F["gcm"];
  VectorXd gcp = F["gcp"];
  VectorXd phi = F["phi"];
  int n = y.size();
  
  MatrixXd Z = design_mat1C(gm, gc, gcm, gcp, x);
  MatrixXd Z1;
  MatrixXd Z2;
  
  if (phi.mean() == 1) {
    Z1 = Z;
    Z2 = Z;
  }
  else {
    Z1 = design_mat2_1C(gm, gc, gcm, gcp, x, phi);
    Z2 = design_mat2_2C(gm, gc, gcm, gcp, x, phi);
  }
  
  double P1 = 1 - pa;
  double P2 = pa;
  
  VectorXd y1 = Z1 * beta;
  VectorXd y2 = Z2 * beta;
  
  VectorXd betagr1(n);
  VectorXd betagr2(n);
  for (int u = 0; u < n; u++) {
    betagr1[u] = P1 * (2*y[u] - 1) * exp(y1[u]) / pow(1 + exp(y1[u]), 2);
    betagr2[u] = P2 * (2*y[u] - 1) * exp(y2[u]) / pow(1 + exp(y2[u]), 2);
  }
  
  MatrixXd betagr(n,6);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < 5; v++) {
      betagr(u,v) = Z1(u,v) * betagr1[u] + Z2(u,v) * betagr2[u];
    }
    betagr(u,5) = exp(y[u]*y2[u]) / (1 + exp(y2[u])) - exp(y[u]*y1[u]) / (1 + exp(y1[u])); 
  }
  return betagr;
}

// [[Rcpp::export]]
MatrixXd log_prY0_grC(VectorXd y, VectorXd gm, VectorXd gc, VectorXd x, VectorXd beta, double pa) {
  List F = distinguish_snpC(gm, gc);
  VectorXd gcm = F["gcm"];
  VectorXd gcp = F["gcp"];
  VectorXd phi = F["phi"];
  int n = y.size();
  
  MatrixXd Z = design_mat1C(gm, gc, gcm, gcp, x);
  MatrixXd Z1;
  MatrixXd Z2;
  
  if (phi.mean() == 1) {
    Z1 = Z;
    Z2 = Z;
  }
  else {
    Z1 = design_mat2_1C(gm, gc, gcm, gcp, x, phi);
    Z2 = design_mat2_2C(gm, gc, gcm, gcp, x, phi);
  }
  
  double P1 = 1 - pa;
  double P2 = pa;
  
  VectorXd y1 = Z1 * beta;
  VectorXd y2 = Z2 * beta;
  
  VectorXd res(n);
  for (int u = 0; u < n; u++) {
    res[u] = P1 * exp(y[u]*y1[u]) / (1 + exp(y1[u])) + P2 * exp(y[u]*y2[u]) / (1 + exp(y2[u]));
  }
  
  VectorXd betagr1(n);
  VectorXd betagr2(n);
  for (int u = 0; u < n; u++) {
    betagr1[u] = P1 * (2*y[u] - 1) * exp(y1[u]) / (pow(1 + exp(y1[u]), 2) * res[u]);
    betagr2[u] = P2 * (2*y[u] - 1) * exp(y2[u]) / (pow(1 + exp(y2[u]), 2) * res[u]);
  }
  
  MatrixXd betagr(n,6);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < 5; v++) {
      betagr(u,v) = Z1(u,v) * betagr1[u] + Z2(u,v) * betagr2[u];
    }
    betagr(u,5) = (exp(y[u]*y2[u]) / (1 + exp(y2[u])) - exp(y[u]*y1[u]) / (1 + exp(y1[u]))) / res[u]; 
  }
  return betagr;
}

// Pr(D=d|gcc,y,x;delta)
// [[Rcpp::export]]
VectorXd prDC(VectorXd d, VectorXd gm, VectorXd gc, VectorXd y, VectorXd x, VectorXd delta) {
  int n = d.size();
  int p = delta.size();
  MatrixXd Z(n,p);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < p; v++) {
      Z(u,0) = 1;
      Z(u,1) = gm[u];
      Z(u,2) = gc[u];
      Z(u,3) = y[u];
      Z(u,4) = x[u];
    }
  }
  VectorXd dbase = Z * delta;
  VectorXd res(n);
  for (int u = 0; u < n; u++) {
    res[u] = exp(d[u] * dbase[u]) / (1 + exp(dbase[u]));
  }
  return res;
}

MatrixXd prD_grC(VectorXd d, VectorXd gm, VectorXd gc, VectorXd y, VectorXd x, VectorXd delta) {
  int n = d.size();
  int p = delta.size();
  MatrixXd Z(n,p);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < p; v++) {
      Z(u,0) = 1;
      Z(u,1) = gm[u];
      Z(u,2) = gc[u];
      Z(u,3) = y[u];
      Z(u,4) = x[u];
    }
  }
  VectorXd dbase = Z * delta;
  MatrixXd resgr(n,p);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < p; v++) {
      resgr(u,v) = Z(u,v) * (2*d[u] - 1) * exp(dbase[u]) / pow(1 + exp(dbase[u]), 2);
    }
  }
  return resgr;
}

MatrixXd log_prD_grC(VectorXd d, VectorXd gm, VectorXd gc, VectorXd y, VectorXd x, VectorXd delta) {
  int n = d.size();
  int p = delta.size();
  MatrixXd Z(n,p);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < p; v++) {
      Z(u,0) = 1;
      Z(u,1) = gm[u];
      Z(u,2) = gc[u];
      Z(u,3) = y[u];
      Z(u,4) = x[u];
    }
  }
  VectorXd dbase = Z * delta;
  
  VectorXd res(n);
  for (int u = 0; u < n; u++) {
    res[u] = exp(d[u] * dbase[u]) / (1 + exp(dbase[u]));
  }
  
  MatrixXd resgr(n,p);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < p; v++) {
      resgr(u,v) = Z(u,v) * (2*d[u] - 1) * exp(dbase[u]) / (pow(1 + exp(dbase[u]), 2) * res[u]);
    }
  }
  return resgr;
}

// H_i=\sum{gc,y}Pr(D=1|gc,y,xi;delta)Pr(Y=y|gm_i,gc,Xi;beta,pa)Pr(gc|gm_i;pa)
double HiC(VectorXd gm, VectorXd x, int i, VectorXd beta, VectorXd delta, double pa) {
  VectorXd gc(6);
  VectorXd y(6);
  VectorXd gm0(6);
  VectorXd x0(6);
  VectorXd d(6);
  gc[0] = 0; gc[1] = 1; gc[2] = 2; gc[3] = 0; gc[4] = 1; gc[5] = 2; 
  y[0] = y[1] = y[2] = 0; y[3] = y[4] = y[5] = 1;
  gm0[0] = gm0[1] = gm0[2] = gm0[3] = gm0[4] = gm0[5] = gm[i-1];
  x0[0] = x0[1] = x0[2] = x0[3] = x0[4] = x0[5] = x[i-1];
  d[0] = d[1] = d[2] = d[3] = d[4] = d[5] = 1;
  
  VectorXd resD = prDC(d, gm0, gc, y, x0, delta);
  VectorXd resY0 = prY0C(y, gm0, gc, x0, beta, pa);
  
  double res = 0;
  MatrixXd mat = gc_on_gmC(pa);
  for (int u = 0; u < 6; u++) {
    res += resD[u] * resY0[u] * mat(gm0[u],gc[u]);
  }
  return res;
}

// d(H_i=\sum{gc,y}Pr(D=1|gc,y,xi;delta)Pr(Y=y|gm_i,gc,Xi;beta,pa)Pr(gc|gm_i;pa))|d(beta,delta,pa)
VectorXd Hi_grC(VectorXd gm, VectorXd x, int i, VectorXd beta, VectorXd delta, double pa) {
  VectorXd gc(6);
  VectorXd y(6);
  VectorXd gm0(6);
  VectorXd x0(6);
  VectorXd d(6);
  gc[0] = 0; gc[1] = 1; gc[2] = 2; gc[3] = 0; gc[4] = 1; gc[5] = 2; 
  y[0] = y[1] = y[2] = 0; y[3] = y[4] = y[5] = 1;
  gm0[0] = gm0[1] = gm0[2] = gm0[3] = gm0[4] = gm0[5] = gm[i-1];
  x0[0] = x0[1] = x0[2] = x0[3] = x0[4] = x0[5] = x[i-1];
  d[0] = d[1] = d[2] = d[3] = d[4] = d[5] = 1;
  
  VectorXd D = prDC(d, gm0, gc, y, x0, delta);
  VectorXd Y0 = prY0C(y, gm0, gc, x0, beta, pa);
  MatrixXd mat = gc_on_gmC(pa);
  MatrixXd Dgr = prD_grC(d, gm0, gc, y, x0, delta);
  MatrixXd Y0gr = prY0_grC(y, gm0, gc, x0, beta, pa);
  MatrixXd matgr = gc_on_gm_grC(pa);
  
  int p = beta.size() + delta.size() + 1;
  VectorXd res(p);
  res.fill(0);
  
  for (int v = 0; v < beta.size(); v++) {
    for (int u = 0; u < 6; u++) {
      res[v] += Y0gr(u,v) * D[u] * mat(gm0[u],gc[u]);
    }
  }
  
  for (int v = 0; v < delta.size(); v++) {
    for (int u = 0; u < 6; u++) {
      res[beta.size()+v] += Dgr(u,v) * Y0[u] * mat(gm0[u],gc[u]);
    }
  }
  
  for (int u = 0; u < 6; u++) {
    res[p-1] += D[u] * (mat(gm0[u],gc[u]) * Y0gr(u,5) + matgr(gm0[u],gc[u]) * Y0[u]);
  }
  
  return res;
}

// functions for single SNP data
// the modified profile likelihood function
// [[Rcpp::export]]
double log_lmpC(VectorXd y, VectorXd d, VectorXd gm, VectorXd gc, VectorXd x, VectorXd beta, VectorXd delta, double pa, double f) {
  int n = y.size();
  int n0 = 0;
  for (int u = 0; u < n; u++) {
    if (d[u] == 0) {
      n0 += 1;
    }
  }
  int n1 = n - n0;
  
  double sum1 = 0;
  VectorXd D = prDC(d, gm, gc, y, x, delta);
  VectorXd Y0 = prY0C(y, gm, gc, x, beta, pa);
  MatrixXd mat = gc_on_gmC(pa);
  for (int u = 0; u < n; u++) {
    sum1 += log(D[u] * Y0[u] * mat(gm[u],gc[u]));
  }
  
  double lambda0 = n1 / (n * f) - n0 / (n * (1 - f));
  double sum2 = 0;
  for (int u = 0; u < n; u++) {
    sum2 += log(1 + lambda0 * (HiC(gm,x,u+1,beta,delta,pa) - f));
  }
  
  double res = sum1 -sum2;
  return res;
}

// gradient of the modified profile likelihood function
// [[Rcpp::export]]
VectorXd log_lmp_grC(VectorXd y, VectorXd d, VectorXd gm, VectorXd gc, VectorXd x, VectorXd beta, VectorXd delta, double pa, double f) {
  int n = y.size();
  int n0 = 0;
  for (int u = 0; u < n; u++) {
    if (d[u] == 0) {
      n0 += 1;
    }
  }
  int n1 = n - n0;
  double lambda0 = n1 / (n * f) - n0 / (n * (1 - f));
  
  MatrixXd betagr = log_prY0_grC(y, gm, gc, x, beta, pa);
  MatrixXd deltagr = log_prD_grC(d, gm, gc, y, x, delta);
  MatrixXd mat = gc_on_gmC(pa);
  MatrixXd matgr = gc_on_gm_grC(pa);
  VectorXd pagr(n);
  for (int u = 0; u < n; u++) {
    pagr[u] = betagr(u,5) + matgr(gm[u],gc[u]) / mat(gm[u],gc[u]);
  }
   
  int p = beta.size() + delta.size() + 1;
  MatrixXd hgr(n,p);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < p; v++) {
      hgr(u,v) = lambda0 * Hi_grC(gm,x,u+1,beta,delta,pa)[v] / (1 + lambda0 * (HiC(gm,x,u+1,beta,delta,pa) - f));
    }
  }
  
  MatrixXd res(n,p);
  for (int v = 0; v < beta.size(); v++) {
    for (int u = 0; u < n; u++) {
      res(u,v) = betagr(u,v) - hgr(u,v);
    }
  }
  
  for (int v = 0; v < delta.size(); v++) {
    for (int u = 0; u < n; u++) {
      res(u,beta.size()+v) = deltagr(u,v) - hgr(u,beta.size()+v);
    }
  }
  
  for (int u = 0; u < n; u++) {
    res(u,p-1)= pagr[u] - hgr(u,p-1);
  }
  
  return res.colwise().sum();
}

// covariance of gradient of the modified profile likelihood function
// [[Rcpp::export]]
MatrixXd cov_log_lmp_grC(VectorXd y, VectorXd d, VectorXd gm, VectorXd gc, VectorXd x, VectorXd beta, VectorXd delta, double pa, double f) {
  int n = y.size();
  int n0 = 0;
  for (int u = 0; u < n; u++) {
    if (d[u] == 0) {
      n0 += 1;
    }
  }
  int n1 = n - n0;
  double lambda0 = n1 / (n * f) - n0 / (n * (1 - f));
  
  MatrixXd betagr = log_prY0_grC(y, gm, gc, x, beta, pa);
  MatrixXd deltagr = log_prD_grC(d, gm, gc, y, x, delta);
  MatrixXd mat = gc_on_gmC(pa);
  MatrixXd matgr = gc_on_gm_grC(pa);
  VectorXd pagr(n);
  for (int u = 0; u < n; u++) {
    pagr[u] = betagr(u,5) + matgr(gm[u],gc[u]) / mat(gm[u],gc[u]);
  }
  
  int p = beta.size() + delta.size() + 1;
  MatrixXd hgr(n,p);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < p; v++) {
      hgr(u,v) = lambda0 * Hi_grC(gm,x,u+1,beta,delta,pa)[v] / (1 + lambda0 * (HiC(gm,x,u+1,beta,delta,pa) - f));
    }
  }
  
  MatrixXd res(n,p);
  for (int v = 0; v < beta.size(); v++) {
    for (int u = 0; u < n; u++) {
      res(u,v) = betagr(u,v) - hgr(u,v);
    }
  }
  
  for (int v = 0; v < delta.size(); v++) {
    for (int u = 0; u < n; u++) {
      res(u,beta.size()+v) = deltagr(u,v) - hgr(u,beta.size()+v);
    }
  }
  
  for (int u = 0; u < n; u++) {
    res(u,p-1)= pagr[u] - hgr(u,p-1);
  }
  
  MatrixXd res1(n0,p);
  MatrixXd res2(n1,p);
  
  int i = 0, j = 0;
  for (int u = 0; u < n; u++) {
    if (d[u] == 0) {
      for (int v = 0; v < p; v++) {
        res1(i,v) = res(u,v);
      }
      i += 1;
    }
    else {
      for (int v = 0; v < p; v++) {
        res2(j,v) = res(u,v);
      }
      j += 1;
    }
  }
  

  MatrixXd t1(n0,p);
  MatrixXd t2(n1,p);
  VectorXd res1mean = res1.colwise().sum()/n0;
  VectorXd res2mean = res2.colwise().sum()/n1;
  
  for (int u = 0; u < n0; u++) {
    for (int v = 0; v < p; v++) {
      t1(u,v) = res1(u,v) - res1mean[v];
    }
  }
  for (int u = 0; u < n1; u++) {
    for (int v = 0; v < p; v++) {
      t2(u,v) = res2(u,v) - res2mean[v];
    }
  }
  
  MatrixXd s1 = (t1.transpose() * t1) * n0 / (n0 - 1);
  MatrixXd s2 = (t2.transpose() * t2) * n1 / (n1 - 1);
  
  return s1+s2;
                                           
}

// functions for haplotype data
// the modified profile likelihood function
// [[Rcpp::export]]
double log_lmp_hapC(VectorXd y, VectorXd d, MatrixXd gmm, MatrixXd gcc, VectorXd x, VectorXd beta, VectorXd delta, int m, MatrixXd hap, VectorXd ppi, double pa, double f) {
  int n = y.size();
  int n0 = 0;
  for (int u = 0; u < n; u++) {
    if (d[u] == 0) {
      n0 += 1;
    }
  }
  int n1 = n - n0;
  
  VectorXd gm = gmm.col(m-1);
  VectorXd gc = gcc.col(m-1);
  
  double sum1 = 0;
  VectorXd D = prDC(d, gm, gc, y, x, delta);
  VectorXd Y = prYC(y, gmm, gcc, x, beta, m, hap, ppi);
  MatrixXd mat = gc_on_gmC(pa);
  for (int u = 0; u < n; u++) {
    sum1 += log(D[u] * Y[u] * mat(gm[u],gc[u]));
  }
  
  double lambda0 = n1 / (n * f) - n0 / (n * (1 - f));
  double sum2 = 0;
  for (int u = 0; u < n; u++) {
    sum2 += log(1 + lambda0 * (HiC(gm,x,u+1,beta,delta,pa) - f));
  }
  
  double res = sum1 -sum2;
  return res;
}

// gradient of the modified profile likelihood function
// [[Rcpp::export]]
VectorXd log_lmp_hap_grC(VectorXd y, VectorXd d, MatrixXd gmm, MatrixXd gcc, VectorXd x, VectorXd beta, VectorXd delta, int m, MatrixXd hap, VectorXd ppi, double pa, double f) {
  int n = y.size();
  int n0 = 0;
  for (int u = 0; u < n; u++) {
    if (d[u] == 0) {
      n0 += 1;
    }
  }
  int n1 = n - n0;
  double lambda0 = n1 / (n * f) - n0 / (n * (1 - f));
  
  VectorXd gm = gmm.col(m-1);
  VectorXd gc = gcc.col(m-1);
  
  MatrixXd betagr = log_prY_grC(y, gmm, gcc, x, beta, m, hap, ppi);
  MatrixXd deltagr = log_prD_grC(d, gm, gc, y, x, delta);
  MatrixXd mat = gc_on_gmC(pa);
  MatrixXd matgr = gc_on_gm_grC(pa);
  VectorXd pagr(n);
  for (int u = 0; u < n; u++) {
    pagr[u] = matgr(gm[u],gc[u]) / mat(gm[u],gc[u]);
  }
  
  int p = beta.size() + delta.size() + 1;
  MatrixXd hgr(n,p);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < p; v++) {
      hgr(u,v) = lambda0 * Hi_grC(gm,x,u+1,beta,delta,pa)[v] / (1 + lambda0 * (HiC(gm,x,u+1,beta,delta,pa) - f));
    }
  }
  
  MatrixXd res(n,p);
  for (int v = 0; v < beta.size(); v++) {
    for (int u = 0; u < n; u++) {
      res(u,v) = betagr(u,v) - hgr(u,v);
    }
  }
  
  for (int v = 0; v < delta.size(); v++) {
    for (int u = 0; u < n; u++) {
      res(u,beta.size()+v) = deltagr(u,v) - hgr(u,beta.size()+v);
    }
  }
  
  for (int u = 0; u < n; u++) {
    res(u,p-1)= pagr[u] - hgr(u,p-1);
  }
  
  return res.colwise().sum();
}

// covariance of gradient of the modified profile likelihood function
// [[Rcpp::export]]
MatrixXd cov_log_lmp_hap_grC(VectorXd y, VectorXd d, MatrixXd gmm, MatrixXd gcc, VectorXd x, VectorXd beta, VectorXd delta, int m, MatrixXd hap, VectorXd ppi, double pa, double f) {
  int n = y.size();
  int n0 = 0;
  for (int u = 0; u < n; u++) {
    if (d[u] == 0) {
      n0 += 1;
    }
  }
  int n1 = n - n0;
  double lambda0 = n1 / (n * f) - n0 / (n * (1 - f));
  
  VectorXd gm = gmm.col(m-1);
  VectorXd gc = gcc.col(m-1);
  
  MatrixXd betagr = log_prY_grC(y, gmm, gcc, x, beta, m, hap, ppi);
  MatrixXd deltagr = log_prD_grC(d, gm, gc, y, x, delta);
  MatrixXd mat = gc_on_gmC(pa);
  MatrixXd matgr = gc_on_gm_grC(pa);
  VectorXd pagr(n);
  for (int u = 0; u < n; u++) {
    pagr[u] = matgr(gm[u],gc[u]) / mat(gm[u],gc[u]);
  }
  
  int p = beta.size() + delta.size() + 1;
  MatrixXd hgr(n,p);
  for (int u = 0; u < n; u++) {
    for (int v = 0; v < p; v++) {
      hgr(u,v) = lambda0 * Hi_grC(gm,x,u+1,beta,delta,pa)[v] / (1 + lambda0 * (HiC(gm,x,u+1,beta,delta,pa) - f));
    }
  }
  
  MatrixXd res(n,p);
  for (int v = 0; v < beta.size(); v++) {
    for (int u = 0; u < n; u++) {
      res(u,v) = betagr(u,v) - hgr(u,v);
    }
  }
  
  for (int v = 0; v < delta.size(); v++) {
    for (int u = 0; u < n; u++) {
      res(u,beta.size()+v) = deltagr(u,v) - hgr(u,beta.size()+v);
    }
  }
  
  for (int u = 0; u < n; u++) {
    res(u,p-1)= pagr[u] - hgr(u,p-1);
  }
  
  MatrixXd res1(n0,p);
  MatrixXd res2(n1,p);
  
  int i = 0, j = 0;
  for (int u = 0; u < n; u++) {
    if (d[u] == 0) {
      for (int v = 0; v < p; v++) {
        res1(i,v) = res(u,v);
      }
      i += 1;
    }
    else {
      for (int v = 0; v < p; v++) {
        res2(j,v) = res(u,v);
      }
      j += 1;
    }
  }
  
  
  MatrixXd t1(n0,p);
  MatrixXd t2(n1,p);
  VectorXd res1mean = res1.colwise().sum()/n0;
  VectorXd res2mean = res2.colwise().sum()/n1;
  
  for (int u = 0; u < n0; u++) {
    for (int v = 0; v < p; v++) {
      t1(u,v) = res1(u,v) - res1mean[v];
    }
  }
  for (int u = 0; u < n1; u++) {
    for (int v = 0; v < p; v++) {
      t2(u,v) = res2(u,v) - res2mean[v];
    }
  }
  
  MatrixXd s1 = (t1.transpose() * t1) * n0 / (n0 - 1);
  MatrixXd s2 = (t2.transpose() * t2) * n1 / (n1 - 1);
  
  return s1+s2;
  
}





