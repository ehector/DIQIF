// [[Rcpp::depends(geepack)]]
#include <Rcpp.h>
#include <math.h>

using namespace std;

#include <tntsupp.h>
#include <geese.h>

extern "C"{
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
}

#include <famstr.h>
#include <param.h>
#include <inter.h>
#include <utils.h>
#include <geesubs.h>


/*********************************************************/

double linkfun_logit(double mu) {return log(mu/(1 - mu));}

double linkinv_logit(double eta) {
  double thres = - log(DBL_EPSILON);
  eta = (eta > thres) ? thres : eta;
  eta = (eta < - thres) ? -thres : eta;
  return exp(eta)/(1 + exp(eta));
}

double mu_eta_logit(double eta) {
  double thres = - log(DBL_EPSILON);
  if (fabs(eta) >= thres) return DBL_EPSILON;
  else return exp(eta)/pow(1 + exp(eta), 2);
} 

bool valideta_logit(double eta) {return true;}

//probit
double linkfun_probit(double mu) {return R::qnorm(mu,0,1,1,0);}
double linkinv_probit(double eta) {
  double thres = -R::qnorm(DBL_EPSILON,0,1,1,0);
  eta = min(thres, max(eta, -thres));
  return R::pnorm(eta,0,1,1,0);
}
double mu_eta_probit(double eta) {
  return max(R::dnorm(eta,0,1,0), DBL_EPSILON);
}
bool valideta_probit(double eta) {return true;}

//cloglog
double linkfun_cloglog(double mu) {return log(-log(1 - mu));}
double linkinv_cloglog(double eta) {
  double ans = 1 - exp(- exp(eta));
  ans = min(1 - DBL_EPSILON, ans);
  return max(DBL_EPSILON, ans);
}
double mu_eta_cloglog(double eta) {
  eta = min(eta, 700.0);
  return max(DBL_EPSILON, exp(eta) * exp(-exp(eta)));
}
bool valideta_cloglog(double eta) {return true;}

//ident
double linkfun_ident(double mu) {return mu;}
double linkinv_ident(double eta) {return eta;}
double mu_eta_ident(double eta) {return 1.0;}
bool valideta_ident(double eta) {return true;}

//log
double linkfun_log(double mu) {return log(mu);}
double linkinv_log(double eta) {return max(DBL_EPSILON, exp(eta));}
double mu_eta_log(double eta) {return max(DBL_EPSILON, exp(eta));}
bool valideta_log(double eta) {return true;}

//sqrt
double linkfun_sqrt(double mu) {return sqrt(mu);}
double inkinv_sqrt(double eta) {return eta * eta;}
double mu_eta_sqrt(double eta) {return 2 * eta;}
bool valideta_sqrt(double eta) {return eta > 0;}

//recipsquare
double linkfun_recipsquare(double mu) {return 1 / mu / mu;}
double linkinv_recipsquare(double eta) {return 1 / sqrt(eta);}
double mu_eta_recipsquare(double eta) {return -1 / (2 * pow(eta, 1.5));}
bool valideta_recipsquare(double eta) {return eta > 0;}

//inverse
double linkfun_inverse(double mu) {return 1 / mu;}
double linkinv_inverse(double eta) {return 1 / eta;}
double mu_eta_inverse(double eta) {return -1 / eta / eta;}
bool valideta_inverse(double eta) {return eta != 0;}

//fisherz
double linkfun_fisherz(double mu) {return log(2/(1 - mu) - 1);}
double linkinv_fisherz(double eta) {
  double thres = - log(DBL_EPSILON);
  eta = (eta > thres) ? thres : eta;
  eta = (eta < - thres) ? -thres : eta;
  return 1 - 2 / (exp(eta) + 1);
}
double mu_eta_fisherz(double eta) {
  double thres = - log(DBL_EPSILON);
  if (fabs(eta) >= thres) return DBL_EPSILON;
  return 2 * exp(eta) / pow(1 + exp(eta), 2);
}
bool valideta_fisherz(double eta) {return true;}


//Lin, Wei, Ying
double linkfun_lwyBC2(double mu) {
  return log(sqrt(mu + 1) - 1);
}
double linkinv_lwyBC2(double eta) {
  double foo = max(DBL_EPSILON, exp(eta));
  return pow(1 + foo, 2.0) - 1;
}
double mu_eta_lwyBC2(double eta) {
  double foo = exp(eta);
  return max(DBL_EPSILON, 2 * (1 + foo) * foo);
}

double linkfun_lwylog(double mu) {
  return log(exp(mu) - 1);
}
double linkinv_lwylog(double eta) {
  return log(exp(eta) + 1);
}
double mu_eta_lwylog(double eta) {
  double foo = exp(eta);
  return foo/(foo + 1);
}

//variance functions
double variance_binomial(double mu) {return mu * (1 - mu);}
double v_mu_binomial(double mu) {return 1 - 2 * mu;}
bool validmu_binomial(double mu) {return mu > 0 && mu < 1;}

double variance_gaussian(double mu) {return 1.0;}
double v_mu_gaussian(double mu) {return .0;}
bool validmu_gaussian(double mu) {return true;}

double variance_poisson(double mu) {return mu;}
double v_mu_poisson(double mu) {return 1.0;}
bool validmu_poisson(double mu) {return mu > 0;}

double variance_inverse_gaussian(double mu) {return pow(mu, 3);}
double v_mu_inverse_gaussian(double mu) {return 3 * mu * mu;}
bool validmu_inverse_gaussian(double mu) {return true;}

double variance_Gamma(double mu) {return mu * mu;}
double v_mu_Gamma(double mu) {return 2 * mu;}
bool validmu_Gamma(double mu) {return mu > 0;}

DMatrix cor_exch(const DVector &rho, const DVector &wave) {
  int n = wave.size();
  DMatrix ans(n,n);
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++) 
      ans(i,j) = (i == j) ? 1.0 : rho(1);
  return ans;
}

DMatrix cor_rho_exch(const DVector &rho, const DVector &wave) {
  int n = wave.size();
  DMatrix ans(n * (n - 1) / 2, 1);
  ans = 1.0;
  return ans;
}

DMatrix cor_indep(const DVector &, const DVector &wave) {
  return ident(wave.size());
}

DMatrix cor_rho_indep(const DVector &, const DVector &) {
  return ident(0);
}

DMatrix cor_fixed(const DVector &rho, const DVector &wave) {
  return cor_unstr(rho, wave);
}

DMatrix cor_rho_fixed(const DVector &, const DVector &) {
  return ident(0);
}


DMatrix cor_ar1(const DVector &rho, const DVector &wave) {
  int n = wave.size();
  DMatrix ans(n,n);
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      ans(i,j) = (i == j) ? 1.0 : pow(rho(1), fabs(wave(j) - wave(i)));
  return ans;
}

DMatrix cor_rho_ar1(const DVector &rho, const DVector &wave) {
  int n = wave.size();
  DMatrix ans(n * (n - 1) / 2, 1);
  int k = 1;
  for (int i = 1; i <= n - 1; i++) {
    for (int j = i + 1; j <= n; j ++) {
      double tmp = fabs(wave(j) - wave(i));
      ans(k, 1) = (tmp == 1.0) ? 1.0 : (tmp * pow(rho(1), tmp - 1.0));
      k++;
    }
  }
  return ans;
}

DMatrix cor_unstr(const DVector &rho, const DVector &wave) {
  DMatrix fullmat = rho2mat(rho);
  return MatRowCol(fullmat, wave, wave);
}

DMatrix cor_rho_unstr(const DVector &rho, const DVector &wave) {
  int n = wave.size();
  return ident(n * (n - 1) / 2);
}

//class Corr
Corr:: Corr(int corst, int maxwave): _corst(corst), _maxwave(maxwave) {
  switch(corst) {
  case INDEPENDENCE:
    _nparam = 0; init(cor_indep, cor_rho_indep); break;
  case EXCHANGEABLE:
    _nparam = 1; init(cor_exch, cor_rho_exch); break;
  case AR1:
    _nparam = 1; init(cor_ar1, cor_rho_ar1); break;
  case UNSTRUCTURED:
  case USERDEFINED:
    _nparam = maxwave; init(cor_unstr, cor_rho_unstr); break;
  case FIXED:
    _nparam = 0; init(cor_fixed, cor_rho_fixed); break;
  }
}

//class Link
//Link::Link() { Link(IDENT); }
//Link::Link(int link) {
Link::Link(int link) {
  switch(link) {
  case LOGIT:
    init(linkfun_logit, linkinv_logit, mu_eta_logit);
    break;
  case IDENT:
    init(linkfun_ident, linkinv_ident, mu_eta_ident);
    break;
  case PROBIT:
    init(linkfun_probit, linkinv_probit, mu_eta_probit);
    break;
  case CLOGLOG:
    init(linkfun_cloglog, linkinv_cloglog, mu_eta_cloglog);
    break;
  case LOG:
    init(linkfun_log, linkinv_log, mu_eta_log);
    break;
  case INVERSE:
    init(linkfun_inverse, linkinv_inverse, mu_eta_inverse);
    break;
  case FISHERZ:
    init(linkfun_fisherz, linkinv_fisherz, mu_eta_fisherz);
    break;
  case LWYBC2:
    init(linkfun_lwyBC2, linkinv_lwyBC2, mu_eta_lwyBC2);
    break;
  case LWYLOG:
    init(linkfun_lwylog, linkinv_lwylog, mu_eta_lwylog);
    break;
  }
}

Link::Link(fun1* linkfun, fun1* linkinv, fun1* mu_eta) {
  init(linkfun, linkinv, mu_eta);
}

//class Variance
//Variance::Variance() {Variance(GAUSSIAN); }
Variance::Variance(int var) {
  //Variance::Variance(int var) {
  switch(var) {
  case GAUSSIAN:
    init(variance_gaussian, v_mu_gaussian, validmu_gaussian); break;
  case BINOMIAL:
    init(variance_binomial, v_mu_binomial, validmu_binomial); break;
  case POISSON:
    init(variance_poisson, v_mu_poisson, validmu_poisson); break;
  case GAMMA:
    init(variance_Gamma, v_mu_Gamma, validmu_Gamma); break;
  }
}

//class GeeStr
GeeStr::GeeStr(int n, Vector<int> meanlink, Vector<int> v,
               Vector<int> scalelink, int corrlink, int scalefix) :    
  CorrLink(corrlink), ScaleFix_(scalefix) {
  //int n = meanlink.size();
  //MeanLink.newsize(n); V.newsize(n); ScaleLink.newsize(n);
  Vector<Link> ML(n), SL(n); Vector<Variance> VS(n);
  for (int i = 1; i <= n; i++) {
    Link ml(meanlink(i)), sl(scalelink(i)); Variance vi(v(i));
    ML(i) = ml; //MeanLink(i) = LINK[meanlink(i) - 1];
    VS(i) = vi; //V(i) = VARIANCE[v(i) - 1]; 
    SL(i) = sl; //ScaleLink(i) = LINK[scalelink(i) - 1];
  }
  MeanLink = ML; V = VS; ScaleLink = SL;
}
DVector GeeStr::MeanLinkfun(const DVector &Mu, const IVector &Wave) {
  int size = Mu.size(); DVector ans(size);
  for (int i = 1; i <= size; i++) ans(i) = MeanLink(Wave(i)).linkfun(Mu(i));
  return ans;
}
DVector GeeStr::MeanLinkinv(const DVector &Eta, const IVector &Wave) {
  int size = Eta.size(); DVector ans(size);
  for (int i = 1; i <= size; i++) ans(i) = MeanLink(Wave(i)).linkinv(Eta(i));
  return ans;
}
DVector GeeStr::MeanMu_eta(const DVector &Eta, const IVector &Wave) {
  int size = Eta.size(); DVector ans(size);
  for (int i = 1; i <= size; i++) ans(i) = MeanLink(Wave(i)).mu_eta(Eta(i));
  return ans;
}
DVector GeeStr::ScaleLinkfun(const DVector &Mu, const IVector &Wave) {
  int size = Mu.size(); DVector ans(size);
  for (int i = 1; i <= size; i++) ans(i) = ScaleLink(Wave(i)).linkfun(Mu(i));
  return ans;
}
DVector GeeStr::ScaleLinkinv(const DVector &Eta, const IVector &Wave) {
  int size = Eta.size(); DVector ans(size);
  for (int i = 1; i <= size; i++) ans(i) = ScaleLink(Wave(i)).linkinv(Eta(i));
  return ans;
}
DVector GeeStr::ScaleMu_eta(const DVector &Eta, const IVector &Wave) {
  int size = Eta.size(); DVector ans(size);
  for (int i = 1; i <= size; i++) ans(i) = ScaleLink(Wave(i)).mu_eta(Eta(i));
  return ans;
}
DVector GeeStr::CorrLinkfun(const DVector &Mu) {
  int size = Mu.size(); DVector ans(size);
  for (int i = 1; i <= size; i++) ans(i) = CorrLink.linkfun(Mu(i));
  return ans;
}
DVector GeeStr::CorrLinkinv(const DVector &Eta) {
  int size = Eta.size(); DVector ans(size);
  for (int i = 1; i <= size; i++) ans(i) = CorrLink.linkinv(Eta(i));
  return ans;
}
DVector GeeStr::CorrMu_eta(const DVector &Eta) {
  int size = Eta.size(); DVector ans(size);
  for (int i = 1; i <= size; i++) ans(i) = CorrLink.mu_eta(Eta(i));
  return ans;
}
DVector GeeStr::v(const DVector &Mu, const IVector &Wave) {
  int size = Mu.size(); DVector ans(size);
  for (int i = 1; i <= size; i++) ans(i) = V(Wave(i)).v(Mu(i));
  return ans;
}
DVector GeeStr::v_mu(const DVector &Mu, const IVector &Wave) {
  int size = Mu.size(); DVector ans(size);
  for (int i = 1; i <= size; i++) ans(i) = V(Wave(i)).v_mu(Mu(i));
  return ans;
}
bool GeeStr::validMu(const DVector &Mu, const IVector &Wave) {
  int size = Mu.size(); 
  bool ans = true;
  for (int i = 1; i <= size; i++) {
    if ( !( V(Wave(i)).validmu(Mu(i)) ) ) {
      ans = false;
      break;
    }
  }
  return ans;
}

/*********************************************************/

IVector comp_lev(GeeStr &geestr, Corr &cor) {
  IVector level(2);
  if (geestr.ScaleFix() != 1) level(1) = 1;
  if (cor.nparam() > 0) level(2) = 1;
  return level;
}

DMatrix gee_infls(DVector &Y, DMatrix &X, 
                  DVector &Offset, DVector &Doffset, DVector &W,
                  IVector &LinkWave, 
                  DMatrix &Zsca, DMatrix &Zcor, DVector &CorP, 
                  IVector &Clusz,  
                  GeeStr &geestr, Corr &cor, GeeParam &par, Control &con) {
  Hess Hi(par), H(par); Grad Gi(par);
  
  int n = Clusz.size();
  IVector ZcorSize(n);
  //if (cor.nparam() > 1) 
  if (cor.corst() > AR1) // == UNSTRUCTRUED || USERDEFINED || FIXED
    for (int i = 1; i <= n; i++) 
      ZcorSize(i) = Clusz(i) * (Clusz(i) - 1) / 2;
  else ZcorSize = 1;
  
  IVector level(2); level = 0;
  if (geestr.ScaleFix() != 1) level(1) = 1;
  if (cor.nparam() > 0) level(2) = 1;
  
  int p = par.p(), q = par.q(), r = par.r();
  DMatrix L11(p,p), L12(p,r), L13(p,q), L22(r,r), L23(r,q), L33(q,q); 
  int l = p + q + r;
  DMatrix infls(l, n), HH(l, l);
  
  Index1D I(0,0), J(0,0);
  Index1D I1(0, 0), JJ(0, 0), I2(0, 0), I3(0, 0);
  I1 = Index1D(1, p);
  I2 = Index1D(p + 1, p + r);
  I3 = Index1D(p + r + 1, p + r + q);
  for (int i = 1; i <= n; i++) {
    int s1 = Clusz(i), s2 = ZcorSize(i), crs = s1 * (s1 - 1) / 2;;
    I = Index1D(1, s1) + I.ubound();
    if (s2 > 0) J = Index1D(1, s2) + J.ubound();
    DVector PRi(s1), Vi(s1), V_Mui(s1); DMatrix Di(s1,p);
    gee_prep(Y, X, Offset, I, LinkWave, par, geestr, PRi, Di, Vi, V_Mui);
    DVector Phii(s1); DMatrix D2i(s1, r);
    PhiandD2(I, LinkWave, Doffset, Zsca, par, geestr, Phii, D2i);
    DMatrix R(s1, s1), E(crs, q);
    RandE(Zcor, I, J, CorP, par, geestr, cor, R, E);
    //cout << "i = " << i;
    DVector Wi = asVec(VecSubs(W, I));
    HiandGi(PRi, Phii, Di, R, Vi, V_Mui, D2i, E, Wi, level, Hi, Gi);
    //cout << "Hi = " << Hi; cout << "H = " << H;
    //cout << "Gi = " << Gi;
    H.inc(Hi); 
    JJ = Index1D(i, i);
    infls(I1, JJ) = asColMat(Gi.U1());
    if (level(1) == 1)  infls(I2, JJ) = asColMat(Gi.U2());
    if (level(2) == 1)  infls(I3, JJ) = asColMat(Gi.U3());
  }
  Hess Hinv = inv(H, level);
  I1 = Index1D(1, p);
  HH(I1, I1) = Hinv.A();
  if (level(1) == 1) {
    HH(I2, I1) = Hinv.B();
    HH(I2, I2) = Hinv.C();
  }
  if (level(2) == 1) {
    HH(I3, I1) = Hinv.D();
    HH(I3, I3) = Hinv.F();
    if (level(1) == 1) HH(I3, I2) = Hinv.E();
  }
  infls = HH * infls;
  return infls;
}


void gee_var(DVector &Y, DMatrix &X, 
             DVector &Offset, DVector &Doffset, DVector &W,
             IVector &LinkWave, 
             DMatrix &Zsca, DMatrix &Zcor, DVector &CorP, 
             IVector &Clusz, IVector &ZcorSize, 
             GeeStr &geestr, Corr &cor, GeeParam &par, Control &con) {
  Hess Hi(par), H(par); Grad Gi(par);
  
  IVector level(2); level = 0;
  if (geestr.ScaleFix() != 1) level(1) = 1;
  if (cor.nparam() > 0) level(2) = 1;
  
  int p = par.p(), q = par.q(), r = par.r();
  DMatrix L11(p,p), L12(p,r), L13(p,q), L22(r,r), L23(r,q), L33(q,q); 
  
  Index1D I(0,0), J(0,0);
  for (int i = 1; i <= Clusz.size(); i++) {
    int s1 = Clusz(i), s2 = ZcorSize(i), crs = s1 * (s1 - 1) / 2;;
    I = Index1D(1, s1) + I.ubound();
    if (s2 > 0) J = Index1D(1, s2) + J.ubound();
    DVector PRi(s1), Vi(s1), V_Mui(s1); DMatrix Di(s1,p);
    gee_prep(Y, X, Offset, I, LinkWave, par, geestr, PRi, Di, Vi, V_Mui);
    DVector Phii(s1); DMatrix D2i(s1, r);
    PhiandD2(I, LinkWave, Doffset, Zsca, par, geestr, Phii, D2i);
    DMatrix R(s1, s1), E(crs, q);
    RandE(Zcor, I, J, CorP, par, geestr, cor, R, E);
    //cout << "i = " << i;
    DVector Wi = asVec(VecSubs(W, I));
    HiandGi(PRi, Phii, Di, R, Vi, V_Mui, D2i, E, Wi, level, Hi, Gi);
    //cout << "Hi = " << Hi; cout << "H = " << H;
    //cout << "Gi = " << Gi;
    H.inc(Hi); 
    L11 = L11 + outerprod(Gi.U1());
    if (level(1) == 1) {
      L12 = L12 + outerprod(Gi.U1(), Gi.U2());
      L22 = L22 + outerprod(Gi.U2());
    }
    if (level(2) == 1) {
      L13 = L13 + outerprod(Gi.U1(), Gi.U3());
      L33 = L33 + outerprod(Gi.U3());
      if (level(1) == 1) L23 = L23 + outerprod(Gi.U2(), Gi.U3());
    }
  }
  //Vbeta:
  Hess Hinv = inv(H, level);
  par.set_vbeta_naiv(Hinv.A());
  par.set_vbeta(Hinv.A() * L11 * Hinv.A());
  
  //Vgamma:
  if (level(1) == 1) {
    par.set_vgamma((Hinv.B() * L11 + Hinv.C() * transpose(L12))  
                     * transpose(Hinv.B()) + 
                       (Hinv.B() * L12 + Hinv.C() * L22) * Hinv.C());
  }
  
  //Valpha:
  if (level(2) == 1) {
    par.set_valpha_naiv(Hinv.F());
    par.set_valpha_stab(Hinv.F() * L33 * Hinv.F());
    par.set_valpha((Hinv.D() * L11 + Hinv.E() * transpose(L12) + 
      Hinv.F() * transpose(L13)) * transpose(Hinv.D()) + 
      (Hinv.D() * L12 + Hinv.E() * L22 + 
      Hinv.F() * transpose(L23)) * transpose(Hinv.E()) +
      (Hinv.D() * L13 + Hinv.E() * L23 + 
      Hinv.F() * L33) * Hinv.F());
  }
}


double update_beta(DVector &Y, DMatrix &X, DVector &Offset, 
                   DVector &W, DVector &Phi,
                   IVector &LinkWave, DVector &CorP,
                   DMatrix &Zcor,  IVector &Clusz,
                   IVector &ZcorSize, IVector &Jack, 
                   GeeParam &par, GeeStr &geestr, Corr &cor) {
  double del = 0;
  //  DVector alp = par.alpha();
  int p = par.p(); 
  DMatrix H(p,p); DVector G(p);
  int n = Clusz.size();
  Index1D I(0,0), J(0,0);
  for (int i = 1; i <= n; i++) {
    int s1 = Clusz(i), s2 = ZcorSize(i);
    I = Index1D(1, s1) + I.ubound(); 
    if (s2 > 0) J = Index1D(1, s2) + J.ubound(); //?? what is s2 == 0 ??
    if (Jack(i) == 1) continue;
    DVector PRi(s1); DMatrix Di(s1,p);
    PRandD(Y, X, Offset, I, LinkWave, par, geestr, PRi, Di);
    DVector rootInvPhii = sqrt(recip(asVec(VecSubs(Phi, I))));
    DVector rootWi = sqrt(asVec(VecSubs(W, I)));
    Di = SMult(rootWi, Di); PRi = SMult(rootWi, PRi);
    Di = SMult(rootInvPhii, Di); PRi = SMult(rootInvPhii, PRi);
    DMatrix R = getR(Zcor, I, J, CorP, par, geestr, cor);
    H = H + AtBiC(Di, R, Di);
    G = G + AtBiC(Di, R, PRi);
  }
  DVector Del = solve(H, G);
  DVector Bnew = par.beta() + Del;
  while (1) {
    //    cerr << "in updating beta: " << "Del = " << Del << endl;
    DVector Eta = X * Bnew + Offset;
    DVector Mu = geestr.MeanLinkinv(Eta, LinkWave);
    if (geestr.validMu(Mu, LinkWave)) break;
    Del = 0.5 * Del;
    Bnew = par.beta() + Del;
  }
  par.set_beta(Bnew);
  del = fmax(fabs(Del));
  return del;
}

double update_gamma(DVector &PR, DVector &W, IVector &LinkWave,
                    IVector &Clusz, IVector &Jack,
                    DVector &Doffset, DMatrix &Zsca, 
                    GeeParam &par, GeeStr &geestr) {
  double del = 0;
  int r = par.r(), n = Clusz.size();
  //  double adj = (double) (PR.size()) / (double)(PR.size() - par.p());
  if (geestr.ScaleFix() == 1) return del; 
  DMatrix H(r,r); DVector G(r);
  Index1D I(0,0);
  for (int i = 1; i <= n; i++) {
    int s1 = Clusz(i);
    I = Index1D(1, s1) + I.ubound();
    if (Jack(i) == 1) continue;
    DVector Phii(s1), Si(s1); DMatrix D2i(s1, r);
    gm_prep(PR, I, LinkWave, Doffset, Zsca, par, geestr, Phii, Si, D2i);
    
    //DMatrix V2 = diag(2.0 * Phii); 
    //independence working structure only now, so no inverting below
    DVector WiV2inv = SMult(asVec(VecSubs(W, I)), recip(2.0 * Phii));
    H = H + Transpose_view(D2i) * SMult(WiV2inv, D2i);
    G = G + Transpose_view(D2i) * SMult(WiV2inv, Si - Phii); //adj * Si
    //H = H + AtBiC(D2i, WiV2, D2i);
    //G = G + AtBiC(D2i, WiV2, Si - Phii);
  }
  DVector Del = solve(H, G);
  //cout << "H = " << H << "G = " << G;
  //par.set_gamma((double) N / (double)(N - p) * (par.gamma() + Del));
  par.set_gamma(par.gamma() + Del);
  del = fmax(fabs(Del));
  return del;
}

double update_alpha(DVector &PR, DVector &Phi, DVector &CorP, DVector &W,
                    IVector &Clusz, IVector &ZcorSize, IVector &Jack,
                    DMatrix &Zcor, 
                    GeeParam &par, GeeStr &geestr, Corr &cor) {
  double del = 0;
  int q = par.q(), n = Clusz.size();
  if (cor.nparam() == 0) return del;
  DMatrix H(q,q); DVector G(q);
  Index1D I(0,0), J(0,0);
  for (int i = 1; i <= n; i++) {
    int s1 = Clusz(i), s2 = ZcorSize(i), crs = s1 * (s1 - 1) / 2;
    I = Index1D(1, s1) + I.ubound(); 
    if (s2 > 0) J = Index1D(1, s2) + J.ubound();
    if (Jack(i) == 1) continue;
    if (s1 == 1) continue;
    
    DVector PRi = asVec(VecSubs(PR, I));
    DVector Phii = asVec(VecSubs(Phi, I));
    DVector sPRi = SMult(reciproot(Phii), PRi);
    DVector zi = genzi(sPRi);
    
    DMatrix R(s1, s1), E(crs, q);
    RandE(Zcor, I, J, CorP, par, geestr, cor, R, E);
    DVector rhoi = utri(R);
    DVector Wi = asVec(VecSubs(W, I));
    //DMatrix V3 = diag(genzi(rootWi)); 
    //independence working correlation only now, no need of inverting below
    DVector WiV3inv = genzi(Wi);
    H = H + Transpose_view(E) * SMult(WiV3inv, E);
    G = G + Transpose_view(E) * SMult(WiV3inv, zi - rhoi);
    //H = H + AtBiC(E, V3, E);
    //G = G + AtBiC(E, V3, zi - rhoi);
  }
  DVector Del = solve(H, G);
  par.set_alpha(par.alpha() + Del);
  del = fmax(fabs(Del));
  return del;
}

/*********************************************************
Input:
Y: response vector;
X: covariate matrix for mean structure;
LinkWave: determines which link to apply on each response component;
Weight: weight, to be implemented ... ... ???;
Offset: offset, to be implemented ... ... ???;
Zsca: covariate matrix for scale structure;
Zcor: covariate matrix for correlation structure;
Corp: correlation parameters to feed cor.mat(rho, .), 
 can be distances for spatial correlation; it is now a vector, which can 
 not really handle >=2 spatial correlations; it really should be a matrix 
 which contains the data to feed cor.mat(rho, .); it actually is the same 
 as LinkWave now, but should be more general to contain high dimensional 
 data, such as coordinates in R x R.
Clusz: cluster sizes;
ZcorSize: number of rows in Zcor for each cluster;
geestr: GEE structure, contains links, variances for each wave;
cor: correlation structure;
par: parameter values;
Jack: Jackknife indicator;
con: control parameters: ScaleFix, ajs, j1s, fij, tol, maxiter;
*********************************************************/
void gee_est(DVector &Y, DMatrix &X, 
             DVector &Offset, DVector &Doffset, DVector &W, 
             IVector &LinkWave,
             DMatrix &Zsca, DMatrix &Zcor, DVector &CorP,
             IVector &Clusz, IVector &ZcorSize,
             GeeStr &geestr, Corr &cor, GeeParam &par,
             IVector &Jack, Control &con) {
  DVector Del(3); 
  
  int N = Y.size();
  DVector PR(N), Phi(N);
  
  int iter; double del;
  for (iter = 0; iter < con.maxiter(); iter++) {
    if (con.trace() == 1) {
      //cerr << "iter " << iter << endl;
      //cerr << "beta = " << par.beta() << "gamma = " << par.gamma() << "alpha = " << par.alpha();
      Rprintf("iter = %d\n", iter);
      Rprintf("beta = "); VecPrint(par.beta()); 
      Rprintf("gamma = "); VecPrint(par.gamma());
      Rprintf("alpha = "); VecPrint(par.alpha());
    }
    //updating beta;
    Phi = getPhi(Doffset, Zsca, LinkWave, par, geestr);
    Del(1) = update_beta(Y, X, Offset, W, Phi, LinkWave, CorP, Zcor, Clusz, ZcorSize, Jack, par, geestr, cor);
    
    //updating gamma;
    PR = getPR(Y, X, Offset, LinkWave, par, geestr);
    //cout << "PR = " << PR;
    //PR = (double) (N / (N - p)) * PR; //df adjusting
    Del(2) = update_gamma(PR, W, LinkWave, Clusz, Jack, Doffset, Zsca, par, geestr);
    
    //updating alpha;
    Phi = getPhi(Doffset, Zsca, LinkWave, par, geestr);
    Del(3) = update_alpha(PR, Phi, CorP, W, Clusz, ZcorSize, Jack, Zcor, par, geestr, cor);
    
    del = fmax(Del);
    if (del <= con.tol()) break;
  }
  if (iter == con.maxiter()) par.set_err(1);
}

void getJackVar(Vector<DVector> &beta_i, Vector<DVector> &alpha_i,
                Vector<DVector> &gamma_i, GeeParam &par, 
                int jack) { //jack = 1, 2, 3 for ajs, j1s, fij
  int I = beta_i.size(), p = par.p(), q = par.q(), r = par.r();
  DMatrix vb(p,p), va(q,q), vc(r,r);
  //cout << par.beta();
  for (int i = 1; i <= I; i++) {
    //cout << "i = " << i << " " << beta_i(i);
    vb = vb + outerprod(beta_i(i) - par.beta());
    //can use level as in gee2_var
    va = va + outerprod(alpha_i(i) - par.alpha());
    vc = vc + outerprod(gamma_i(i) - par.gamma());
  }
  
  double f = (double) (I - p - q - r) / I;
  if (jack == 3) {//fij
    par.set_vbeta_fij(f * vb);
    par.set_valpha_fij(f * va);
    par.set_vgamma_fij(f * vc);
  }
  else if (jack == 2) { //j1s
    par.set_vbeta_j1s(f * vb);
    par.set_valpha_j1s(f * va);
    par.set_vgamma_j1s(f * vc);
  }
  else {//ajs
    par.set_vbeta_ajs(f * vb);
    par.set_valpha_ajs(f * va);
    par.set_vgamma_ajs(f * vc);
  }
}

void gee_jack(DVector &Y, DMatrix &Xmat, DVector &Offset, DVector &Doffset,
              DVector &W, IVector &LinkWave, DMatrix &Zsca, DMatrix &Zcor,
              DVector &CorP, IVector &Clusz, IVector &ZcorSize,
              GeeStr &geestr, Corr &cor,
              GeeParam &par, Control &con) {
  int I = Clusz.size();
  //  int p = par.p(), q = par.q(), r = par.r();
  IVector Jack(I);
  Vector<DVector> beta_i(I), alpha_i(I), gamma_i(I);
  Vector<DVector> beta_fi(I), alpha_fi(I), gamma_fi(I);
  //DVector b0(p), a0(q), c0(r);
  //beta_i = b0; alpha_i(I) = a0; gamma_i(I) = c0;
  //beta_fi = b0; alpha_fi(I) = a0; gamma_fi(I) = c0;
  Control con1(con); con1.set_maxiter(1); //for j1s
  
  for (int i = 1; i <= I; i++) {
    Jack(i) = 1;
    GeeParam par_i(par.beta(), par.alpha(), par.gamma());
    if (con.j1s() == 1) {
      gee_est(Y, Xmat, Offset, Doffset, W, LinkWave, Zsca, Zcor, CorP, Clusz, ZcorSize, geestr, cor, par_i, Jack, con1); //1-step
      beta_i(i) = par_i.beta();
      alpha_i(i) = par_i.alpha();
      gamma_i(i) = par_i.gamma();
    }
    
    if (con.fij() == 1) {
      gee_est(Y, Xmat, Offset, Doffset, W, LinkWave, Zsca, Zcor, CorP, Clusz, ZcorSize, geestr, cor, par_i, Jack, con); //full iterated
      beta_fi(i) = par_i.beta();
      alpha_fi(i) = par_i.alpha();
      gamma_fi(i) = par_i.gamma();
    }
    Jack(i) = 0;
  }
  if (con.j1s() == 1) getJackVar(beta_i, alpha_i, gamma_i, par, 2);
  if (con.fij() == 1) getJackVar(beta_fi, alpha_fi, gamma_fi, par, 3);
}

void jack_ajs(DVector &Y, DMatrix &X, DVector &Offset, DVector &Doffset,
              DVector &W, IVector &LinkWave, DMatrix &Zsca, DMatrix &Zcor,
              DVector &CorP, IVector &Clusz, IVector &ZcorSize,
              GeeStr &geestr, Corr &cor,
              GeeParam &par, Control &con) {//con is not used
  int I = Clusz.size(), p = par.p(), q = par.q(), r = par.r();
  Vector<Hess> His(I); Vector<Grad> Gis(I);
  IVector level = comp_lev(geestr, cor), Scur(Y.size()); Scur = 1;
  HisandGis(Y, X, Offset, Doffset, W, LinkWave, Clusz, ZcorSize,
            Zsca, Zcor, CorP, par, geestr, cor, Scur, level,
            His, Gis);
  Hess Hn(par); 
  for (int i = 1; i <= I; i++) Hn.inc(His(i));
  
  Vector<DVector> beta_i(I), alpha_i(I), gamma_i(I);
  DVector b0(p), a0(q), c0(r);
  beta_i = b0; alpha_i(I) = a0; gamma_i(I) = c0;
  
  DMatrix vb(p,p), va(q,q), vc(r,r);
  for (int i = 1; i <= I; i++) {
    Hess H_i = Hn - His(i);
    H_i = inv(H_i, level);
    beta_i(i) = H_i.A() * Gis(i).U1();
    gamma_i(i) = H_i.B() * Gis(i).U1() + H_i.C() * Gis(i).U2();
    alpha_i(i) = H_i.D() * Gis(i).U1() + H_i.E() * Gis(i).U2() + H_i.F() * Gis(i).U3();
    vb = vb + outerprod(beta_i(i));
    //can use level as in gee2_var
    va = va + outerprod(alpha_i(i));
    vc = vc + outerprod(gamma_i(i));
  }
  double f = (double) (I - p - q - r) / I;
  par.set_vbeta_ajs(f * vb);
  par.set_valpha_ajs(f * va);
  par.set_vgamma_ajs(f * vc);
}

void gee_top(DVector &Y, DMatrix &Xmat, 
             DVector &Offset, DVector &Doffset, DVector &W, 
             IVector &LinkWave,
             DMatrix &Zsca, DMatrix &Zcor, DVector &CorP,
             IVector &Clusz, 
             GeeStr &geestr, Corr &cor, GeeParam &par,
             Control &con) {
  int I = Clusz.size();
  IVector Jack(I), ZcorSize(I);
  //initializing ZcorSize
  //if (cor.nparam() > 1) 
  if (cor.corst() > AR1) // == UNSTRUCTRUED || USERDEFINED || FIXED
    for (int i = 1; i <= I; i++) 
      ZcorSize(i) = Clusz(i) * (Clusz(i) - 1) / 2;
  else ZcorSize = 1;
  
  gee_est(Y, Xmat, Offset, Doffset, W, LinkWave, Zsca, Zcor, CorP, Clusz, ZcorSize, geestr, cor, par, Jack, con);
  
  gee_var(Y, Xmat, Offset, Doffset, W, LinkWave, Zsca, Zcor, CorP, Clusz, ZcorSize, geestr, cor, par, con);
  
  if (con.ajs() == 1) 
    jack_ajs(Y, Xmat, Offset, Doffset, W, LinkWave, Zsca, Zcor, CorP, Clusz, ZcorSize, geestr, cor, par, con);
  
  if (con.j1s() + con.fij() > 0) 
    gee_jack(Y, Xmat, Offset, Doffset, W, LinkWave, Zsca, Zcor, CorP, Clusz, ZcorSize, geestr, cor, par, con);
}

extern "C" {
  SEXP gee_rap(SEXP y, SEXP x, SEXP offset, SEXP doffset, SEXP w,
               SEXP linkwave, SEXP zsca, SEXP zcor, SEXP corp,
               SEXP clusz, SEXP geestr, SEXP cor, SEXP par, SEXP con) {
    DVector Y = asDVector(y), Offset = asDVector(offset), Doffset = asDVector(doffset), W = asDVector(w);
    IVector LinkWave = asIVector(linkwave); 
    DVector CorP = asDVector(corp); 
    DMatrix X = asDMatrix(x), Zsca = asDMatrix(zsca), Zcor = asDMatrix(zcor);
    IVector Clusz = asIVector(clusz); // ZcorSize = asIVector(zcorsize);
    Control Con = asControl(con);   
    GeeParam Par = asGeeParam(par);   
    GeeStr Geestr = asGeeStr(geestr);   
    Corr Cor = asCorr(cor);   
    
    gee_top(Y, X, Offset, Doffset, W, LinkWave, Zsca, Zcor, CorP, Clusz, Geestr, Cor, Par, Con);
    SEXP ans = asSEXP(Par);
    return ans;
  }
  
  /* return the influence functions for parameters */
  SEXP infls_rap(SEXP y, SEXP x, SEXP offset, SEXP doffset, SEXP w,
                 SEXP linkwave, SEXP zsca, SEXP zcor, SEXP corp,
                 SEXP clusz, SEXP geestr, SEXP cor, SEXP par, SEXP con) {
    DVector Y = asDVector(y), Offset = asDVector(offset), Doffset = asDVector(doffset), W = asDVector(w);
    IVector LinkWave = asIVector(linkwave); 
    DVector CorP = asDVector(corp); 
    DMatrix X = asDMatrix(x), Zsca = asDMatrix(zsca), Zcor = asDMatrix(zcor);
    IVector Clusz = asIVector(clusz); // ZcorSize = asIVector(zcorsize);
    Control Con = asControl(con);   
    GeeParam Par = asGeeParam(par);   
    GeeStr Geestr = asGeeStr(geestr);   
    Corr Cor = asCorr(cor);   
    
    DMatrix infls = gee_infls(Y, X, Offset, Doffset, W, LinkWave, Zsca, Zcor, CorP, Clusz, Geestr, Cor, Par, Con);
    SEXP ans = asSEXP(infls);
    return ans;
  }
}

/*********************************************************/

Grad & Grad::operator=(const Grad &G) {
  U1_ = G.U1_; U2_ = G.U2_; U3_ = G.U3_;
  return *this;
}

ostream& operator<<(ostream& s, const Grad &G) {
  s << "U1 = " << G.U1() << "U2 = " << G.U2() << "U3 = " << G.U3(); 
  return s;
}

Hess operator-(Hess &H1, Hess &H2) {
  Hess ans(H1);
  ans.dec(H2);
  return ans;
}

Hess inv(Hess &H, IVector &level) {
  Hess ans(H);
  ans.set_A(inv(H.A()));
  if (level(1) == 1) {
    ans.set_C(inv(H.C()));
    ans.set_B(-1.0 * ans.C() * H.B() * ans.A());
  }
  if (level(2) == 1) {
    ans.set_F(inv(H.F()));
    ans.set_E(-1.0 * ans.F() * H.E() * ans.C());
    ans.set_D(-1.0 * ans.F() * (H.D() * ans.A() + H.E() * ans.B()));
  }
  return ans;
}

Hess operator*(const double &x, const Hess &H) {
  Hess ans(H);
  ans.set_A(x * H.A()); ans.set_B(x * H.B()); ans.set_C(x * H.C());
  ans.set_D(x * H.D()); ans.set_E(x * H.E()); ans.set_F(x * H.F());
  return ans;
}

ostream& operator<<(ostream& s, const Hess &H) {
  s << "A = " << H.A() << "B = " << H.B() << "C = " << H.C() 
    << "D = " << H.D()<< "E = " << H.E() << "F = " << H.F();
  return s;
}

DVector genzi(const DVector &PR) {
  int n = PR.size();
  DVector ans(n * (n - 1)/2);
  int k = 1;
  for (int i = 1; i <= n - 1; i++)
    for (int j = i + 1; j <= n; j++) ans(k++) = PR(i) * PR(j);
  return ans;
}

DVector utri(const DMatrix &R) {
  int n = R.dim(1);
  //assert (n > 1);
  DVector ans(n * (n - 1) / 2);
  int k = 1;
  for (int i = 1; i <= n - 1; i++)
    for (int j = i + 1; j <= n; j++) ans(k++) = R(i,j);
  return ans;
}

DMatrix getZ_Beta(DMatrix &D, DVector &PR, 
                  DVector &V, DVector &V_Mu, DVector &z) {
  //note: this is the version which excludes phi in the formula
  DMatrix ans(z.size(), D.dim(2));
  int k = 1, n = PR.size();
  for (int i = 1; i <= n - 1; i++) {
    DMatrix Di = asMat(MatRow(D,i));
    for (int j = i + 1; j <= n; j++) {
      DMatrix Dj = asMat(MatRow(D,j));
      DMatrix foo = V_Mu(i) * reciproot(V(i)) * Di + 
        V_Mu(j) * reciproot(V(j)) * Dj;
      DMatrix bar = - PR(i) * Di - PR(j) * Dj - 0.5 * PR(i) * PR(j) * foo;
      //cout << "bar = " << bar << "k = " << k;
      MatRow(ans, k) = bar;
      //cout << " ans = " << ans;
      k++;
    }
  }
  return ans;
}

DMatrix getZ_Gamma(DMatrix &D, DVector &PR, DVector &Phi, DVector &z) {
  DMatrix ans(z.size(), D.dim(2));
  int k = 1, n = PR.size();
  for (int i = 1; i <= n - 1; i++) {
    DMatrix Di = asMat(MatRow(D,i));
    for (int j = i + 1; j <= n; j++) {
      DMatrix Dj = asMat(MatRow(D,j));
      //MatRow(ans, k) = -0.5 * z(k) * (sqrt(Phi(j)/Phi(i)) * Di +
      //sqrt(Phi(i)/Phi(j)) * Dj);
      //This has caused the scale problem; The first time, scale problem was caused by operator * for Hess, where one component did not get touched, in the old geese (LAPACK);
      MatRow(ans, k) = -0.5 * z(k) * (1.0 / Phi(i) * Di +
        1.0 / Phi(j) * Dj);
      k++;
    }
  }
  return ans;
}

DMatrix getS_Beta(DMatrix &D, DVector &PR, DVector &V, DVector &V_Mu) {
  DMatrix ans(D);
  for (int i = 1; i <= ans.dim(1); i++) {
    DMatrix Di = asMat(MatRow(D,i));
    double f = -2 * PR(i) / sqrt(V(i)) - PR(i) * PR(i)/V(i) * V_Mu(i);
    MatRow(ans, i) = f * Di;
  }
  return ans;
}

void HiandGi(DVector &PRi, DVector &Phii, DMatrix &Di, DMatrix &R,
             DVector &Vi, DVector &V_Mui, DMatrix &D2i, DMatrix &E,
             DVector &Wi, IVector &level,
             Hess &H, Grad &G) {
  int s = PRi.size();
  //beta
  DVector rootPhii = sqrt(Phii);
  DMatrix V1 = diag(rootPhii) * R * diag(rootPhii);
  DVector rootWi = sqrt(Wi);
  DMatrix rootWD = SMult(rootWi, Di); 
  DVector rootWPR = SMult(rootWi, PRi);
  H.set_A(AtBiC(rootWD, V1, rootWD));
  G.set_U1(AtBiC(rootWD, V1, rootWPR));
  //H.set_A(AtBiC(Di, V1, Di));
  //G.set_U1(AtBiC(Di, V1, PRi));
  //gamma
  if (level(1) == 1) {//if (par.ScaleFix() != 1) {
    DVector Si = square(PRi);
    DVector WiV2inv = SMult(Wi, recip(2.0 * Phii));
    H.set_C(Transpose_view(D2i) * SMult(WiV2inv, D2i));
    DMatrix S_Beta = getS_Beta(Di, PRi, Vi, V_Mui);
    H.set_B(Transpose_view(D2i) * SMult(-1.0 * WiV2inv, S_Beta));
    G.set_U2(Transpose_view(D2i) * SMult(WiV2inv, Si - Phii));
    //DMatrix V2 = diag(2.0 * Phii);
    //H.set_C(AtBiC(D2i, V2, D2i));
    //DMatrix S_Beta = getS_Beta(Di, PRi, Vi, V_Mui);
    //H.set_B(AtBiC(D2i, V2, S_Beta));
    //G.set_U2(AtBiC(D2i, V2, S - Phii));
  }
  //alpha
  if (level(2) == 1) {//if (cor.nparam() > 0) {
    if (s == 1) return;
    DVector sPRi = SMult(reciproot(Phii), PRi);
    DVector zi = genzi(sPRi);
    DVector rhoi = utri(R); 
    //DMatrix W = ident(s * (s - 1) / 2);
    DVector Sca = genzi(reciproot(Phii)); 
    DVector WiV3inv = genzi(Wi);
    
    //H.set_F(AtBiC(E, W, E));
    H.set_F(Transpose_view(E) * SMult(WiV3inv, E));
    DMatrix Z_Beta = getZ_Beta(Di, PRi, Vi, V_Mui, zi);
    Z_Beta = SMult(Sca, Z_Beta);
    //H.set_D(AtBiC(E, W, Z_Beta));
    H.set_D(Transpose_view(E) * SMult(-1.0 * WiV3inv, Z_Beta));
    //G.set_U3(AtBiC(E, W, zi - rhoi));
    G.set_U3(Transpose_view(E) * SMult(WiV3inv, zi - rhoi));
    if (level(1) == 1) {//if (par.ScaleFix() != 1) {
      DMatrix Z_Gamma = getZ_Gamma(D2i, PRi, Phii, zi);
      //H.set_E(AtBiC(E, W, Z_Gamma));
      H.set_E(Transpose_view(E) * SMult(-1.0 * WiV3inv, Z_Gamma));
    }
  }
}

void PRandD(DVector &Yi, DMatrix &Xi, DVector &Offseti,
            IVector &Wavei, GeeParam &par, GeeStr &geestr,
            DVector &PRi, DMatrix &Di) {
  DVector Eta = Xi * par.beta() + Offseti;
  DVector Mu = geestr.MeanLinkinv(Eta, Wavei);
  DVector V = geestr.v(Mu, Wavei);
  DVector Mu_Eta = geestr.MeanMu_eta(Eta, Wavei);
  
  DVector InvRootV = reciproot(V);
  Di = SMult(InvRootV, SMult(Mu_Eta, Xi));
  PRi = SMult(InvRootV, Yi - Mu);
}

void PRandD(DVector &Y, DMatrix &X, DVector &Offset,
            Index1D &I, IVector &LinkWave, 
            GeeParam &par, GeeStr &geestr,
            DVector &PRi, DMatrix &Di) {
  DVector Yi = asVec(VecSubs(Y, I));
  DMatrix Xi = asMat(MatRows(X, I));
  DVector Offseti = asVec(VecSubs(Offset, I));
  IVector Wavei = asVec(VecSubs(LinkWave, I));
  DVector Eta = Xi * par.beta() + Offseti;
  DVector Mu = geestr.MeanLinkinv(Eta, Wavei);
  DVector V = geestr.v(Mu, Wavei);
  DVector Mu_Eta = geestr.MeanMu_eta(Eta, Wavei);
  
  DVector InvRootV = reciproot(V);
  Di = SMult(InvRootV, SMult(Mu_Eta, Xi));
  PRi = SMult(InvRootV, Yi - Mu);
}

void gee_prep(DVector &Yi, DMatrix &Xi, DVector &Offseti,
              IVector &Wavei, GeeParam &par, GeeStr &geestr,
              DVector &PRi, DMatrix &Di, DVector &Vi, DVector &V_Mui) {
  DVector Eta = Xi * par.beta() + Offseti;
  DVector Mu = geestr.MeanLinkinv(Eta, Wavei);
  DVector V = geestr.v(Mu, Wavei);
  DVector Mu_Eta = geestr.MeanMu_eta(Eta, Wavei);
  
  DVector InvRootV = reciproot(V);
  Di = SMult(InvRootV, SMult(Mu_Eta, Xi));
  PRi = SMult(InvRootV, Yi - Mu);
  Vi = geestr.v(Mu, Wavei);
  V_Mui = geestr.v_mu(Mu, Wavei);
}

void gee_prep(DVector &Y, DMatrix &X, DVector &Offset,
              Index1D &I, IVector &LinkWave, 
              GeeParam &par, GeeStr &geestr,
              DVector &PRi, DMatrix &Di, DVector &Vi, DVector &V_Mui) {
  DVector Yi = asVec(VecSubs(Y, I));
  DMatrix Xi = asMat(MatRows(X, I));
  DVector Offseti = asVec(VecSubs(Offset, I));
  IVector Wavei = asVec(VecSubs(LinkWave, I));
  DVector Eta = Xi * par.beta() + Offseti;
  DVector Mu = geestr.MeanLinkinv(Eta, Wavei);
  DVector V = geestr.v(Mu, Wavei);
  DVector Mu_Eta = geestr.MeanMu_eta(Eta, Wavei);
  
  DVector InvRootV = reciproot(V);
  Di = SMult(InvRootV, SMult(Mu_Eta, Xi));
  PRi = SMult(InvRootV, Yi - Mu);
  Vi = geestr.v(Mu, Wavei);
  V_Mui = geestr.v_mu(Mu, Wavei);
}

DMatrix getR(DMatrix &Zmati, DVector &corp, 
             GeeParam &par, GeeStr &geestr, Corr &cor) {
  DVector alp = par.alpha(); 
  int s = corp.dim(); // corp should determine meta par for R
  if (s == 1) return ident(1); 
  else if (cor.corst() == INDEPENDENCE) //indenpendence
    return cor.mat(alp, corp); 
  else {
    DVector Eta = Zmati * alp;
    DVector Rho = geestr.CorrLinkinv(Eta);
    return cor.mat(Rho, corp);
  }
}

DMatrix getR(DMatrix &Zmat, Index1D &I, Index1D &J, DVector &CorP,
             GeeParam &par, GeeStr &geestr, Corr &cor) {
  DVector alp = par.alpha(); 
  DVector corp = asVec(VecSubs(CorP, I)); 
  int s = corp.dim(); // corp should determine meta par for R
  if (s == 1) return ident(1); 
  else if (cor.corst() == INDEPENDENCE) //indenpendence
    return cor.mat(alp, corp); 
  else{
    DMatrix Zmati = asMat(MatRows(Zmat, J));
    DVector Eta = Zmati * alp;
    DVector Rho = geestr.CorrLinkinv(Eta);
    return cor.mat(Rho, corp);
  }
}

int RandE(DMatrix &Zmati, DVector &corp,
          GeeParam &par, GeeStr &geestr, Corr &cor,
          DMatrix &R, DMatrix &E) {
  DVector alp = par.alpha();
  //DVector corp = asVec(VecSubs(CorP, I));
  int s = corp.dim();
  if (s == 1) {
    R = ident(1); 
    return 0;
  }
  else if (cor.corst() == INDEPENDENCE) { //no need for E
    R = cor.mat(alp, corp);
    return 0;
  }
  else if (cor.corst() == FIXED) {
    DVector Eta = Zmati * alp;
    DVector Rho = geestr.CorrLinkinv(Eta);
    R = cor.mat(Rho, corp);
    return 0;
  }
  else {
    //DMatrix Zmati = asMat(MatRows(Zmat, J)); 
    DVector Eta = Zmati * alp;
    DVector Rho = geestr.CorrLinkinv(Eta);
    R = cor.mat(Rho, corp);
    DVector Rho_Alp = geestr.CorrMu_eta(Eta);
    DMatrix Cor_Rho = cor.cor_rho(Rho, corp);
    E = Cor_Rho * SMult(Rho_Alp,  Zmati);
    return 0;
  }
}

int RandE(DMatrix &Zmat, Index1D &I, Index1D &J, DVector &CorP,
          GeeParam &par, GeeStr &geestr, Corr &cor,
          DMatrix &R, DMatrix &E) {
  DVector alp = par.alpha();
  DVector corp = asVec(VecSubs(CorP, I));
  int s = corp.dim();
  if (s == 1) {
    R = ident(1); 
    return 0;
  }
  else if (cor.corst() == INDEPENDENCE) { //no need for E
    R = cor.mat(alp, corp);
    return 0;
  }
  else if (cor.corst() == FIXED) {
    DMatrix Zmati = asMat(MatRows(Zmat, J)); 
    DVector Eta = Zmati * alp;
    DVector Rho = geestr.CorrLinkinv(Eta);
    R = cor.mat(Rho, corp);
    return 0;
  }
  else {
    DMatrix Zmati = asMat(MatRows(Zmat, J)); 
    DVector Eta = Zmati * alp;
    DVector Rho = geestr.CorrLinkinv(Eta);
    R = cor.mat(Rho, corp);
    DVector Rho_Alp = geestr.CorrMu_eta(Eta);
    DMatrix Cor_Rho = cor.cor_rho(Rho, corp);
    E = Cor_Rho * SMult(Rho_Alp,  Zmati);
    return 0;
  }
}

void gm_prep(DVector &PRi, IVector &Wavei,
             DVector &Doffseti, DMatrix &Zi, GeeParam &par, GeeStr &geestr,
             DVector &Phii, DVector &Si, DMatrix &D2i) {
  DVector Zeta = Zi * par.gamma() + Doffseti;
  DVector Phi_Zeta = geestr.ScaleMu_eta(Zeta, Wavei);
  
  Phii = geestr.ScaleLinkinv(Zeta, Wavei);
  Si = square(PRi);
  D2i = Phi_Zeta * Zi;
}

void gm_prep(DVector &PR, Index1D &I, IVector &LinkWave,
             DVector &Doffset, DMatrix &Zsca, GeeParam &par, GeeStr &geestr,
             DVector &Phii, DVector &Si, DMatrix &D2i) {
  DMatrix Zi = asMat(MatRows(Zsca, I));
  DVector Doffseti = asVec(VecSubs(Doffset, I));
  IVector Wavei = asVec(VecSubs(LinkWave, I));
  DVector Zeta = Zi * par.gamma() + Doffseti;
  DVector Phi_Zeta = geestr.ScaleMu_eta(Zeta, Wavei);
  DVector PRi = asVec(VecSubs(PR, I));
  
  Phii = geestr.ScaleLinkinv(Zeta, Wavei);
  Si = square(PRi);
  D2i = Phi_Zeta * Zi;
}

void PhiandD2(IVector &Wavei,
              DVector &Doffseti, DMatrix &Zi, GeeParam &par, GeeStr &geestr,
              DVector &Phii, DMatrix &D2i) {
  DVector Zeta = Zi * par.gamma() + Doffseti;
  Phii = geestr.ScaleLinkinv(Zeta, Wavei);
  if (geestr.ScaleFix() == 1) return;
  DVector Phi_Zeta = geestr.ScaleMu_eta(Zeta, Wavei);
  D2i = Phi_Zeta * Zi;
}

void PhiandD2(Index1D &I, IVector &LinkWave,
              DVector &Doffset, DMatrix &Zsca, GeeParam &par, GeeStr &geestr,
              DVector &Phii, DMatrix &D2i) {
  DMatrix Zi = asMat(MatRows(Zsca, I));
  DVector Doffseti = asVec(VecSubs(Doffset, I));
  IVector Wavei = asVec(VecSubs(LinkWave, I));
  DVector Zeta = Zi * par.gamma() + Doffseti;
  Phii = geestr.ScaleLinkinv(Zeta, Wavei);
  if (geestr.ScaleFix() == 1) return;
  DVector Phi_Zeta = geestr.ScaleMu_eta(Zeta, Wavei);
  D2i = Phi_Zeta * Zi;
}

DVector getPR(DVector &Y, DMatrix &X, DVector &Offset, IVector &LinkWave,
              GeeParam &par, GeeStr &geestr) {
  DVector Eta = X * par.beta() + Offset;
  DVector Mu = geestr.MeanLinkinv(Eta, LinkWave);
  DVector V = geestr.v(Mu, LinkWave);
  DVector InvRootV = reciproot(V);
  return SMult(InvRootV, Y - Mu);
}

DVector getPhi(DVector &Doffset, DMatrix &Zsca, IVector &LinkWave,
               GeeParam &par, GeeStr &geestr) {
  DVector Zeta = Zsca * par.gamma() + Doffset;
  return geestr.ScaleLinkinv(Zeta, LinkWave);
}


void getDatI(DVector &Y, DVector &Offset, DVector &Doffset, 
             DVector &W, DVector &CorP, 
             DMatrix &X, DMatrix &Zsca, DMatrix &Zcor,
             IVector &LinkWave, 
             //extract indicator
             Index1D &I, Index1D &J, IVector Scuri,
             Corr &cor,
             //output
             DVector &VYi, DVector &VOffseti, DVector &VDoffseti, 
             DVector &VWi, DVector &VCorPi,
             DMatrix &VXi, DMatrix &VZscai, DMatrix &VZcori,
             IVector &VLinkWavei) {
  int s = Scuri.size();
  //get dat i
  DVector Yi = asVec(VecSubs(Y, I));
  DVector Offseti = asVec(VecSubs(Offset, I));
  DVector Wi = asVec(VecSubs(W, I));
  DVector CorPi = asVec(VecSubs(CorP, I));
  DMatrix Xi = asMat(MatRows(X, I));
  DMatrix Zscai = asMat(MatRows(Zsca, I));
  IVector LinkWavei = asVec(VecSubs(LinkWave, I));
  DMatrix Zcori; DVector Doffseti;
  if (cor.corst() > INDEPENDENCE && s > 1 )  {
    Zcori = asMat(MatRows(Zcor, J));
  }
  Doffseti = asVec(VecSubs(Doffset, I));
  
  //valid dat i
  IVector VI = genVI(Scuri), VJ = genCrossVI(Scuri);
  VYi = Valid(Yi, VI); VOffseti = Valid(Offseti, VI); 
  VWi = Valid(Wi, VI); VCorPi = Valid(CorPi, VI);
  VXi = Valid(Xi, VI);
  VZscai = Valid(Zscai, VI);
  VLinkWavei = Valid(LinkWavei, VI);
  if (cor.corst() > INDEPENDENCE && s > 1) {
    if (cor.nparam() == 1) VZcori = Zcori;
    else VZcori = Valid(Zcori, VJ);
    //VDoffseti = Valid(Doffseti, VJ); //this is for log odds for ordinal
  }
  VDoffseti = Valid(Doffseti, VI);
}

void HnandGis(DVector &Y, DMatrix &X, 
              DVector &Offset, DVector &Doffset, DVector &W, 
              IVector &LinkWave, IVector &Clusz, IVector &ZcorSize,
              DMatrix &Zsca, DMatrix &Zcor, DVector &CorP,
              GeeParam &par, GeeStr &geestr, Corr &cor,
              IVector &Scur, IVector &level, //Scur is the valid data indicator
              //output
              Hess &Hn, Vector<Grad> &Gis) {
  int N = Clusz.size();
  Hess H(par);
  Vector<Hess> His(N); His = H; 
  HisandGis(Y, X, Offset, Doffset, W, LinkWave, Clusz, ZcorSize,
            Zsca, Zcor, CorP, par, geestr, cor, Scur, level,
            His, Gis);
  for (int i = 1; i <= N; i++) H.inc(His(i));
  Hn = (1.0/(double) N) * H;
}

void HisandGis(DVector &Y, DMatrix &X, 
               DVector &Offset, DVector &Doffset, DVector &W, 
               IVector &LinkWave, IVector &Clusz, IVector &ZcorSize,
               DMatrix &Zsca, DMatrix &Zcor, DVector &CorP,
               GeeParam &par, GeeStr &geestr, Corr &cor,
               IVector &Scur, IVector &level,
               //output
               Vector<Hess> &His, Vector<Grad> &Gis) {
  Index1D I(0,0), J(0,0);
  int N = Clusz.size();
  int pb = par.p(), pg = par.r(), pa = par.q();
  DVector V0(pb);
  Hess H(par), Hi(par); Grad Gi(par);
  //cout << "N = " << N;
  for (int i = 1; i <= N; i++) {
    int s1 = Clusz(i), s2 = ZcorSize(i);
    I = Index1D(1, s1) + I.ubound();
    if (s2 > 0) J = Index1D(1, s2) + J.ubound();
    
    IVector Scuri = asVec(VecSubs(Scur, I));
    if (sum(Scuri) == 0)  continue;
    //get and valid data i
    DVector Yi, Offseti, Doffseti, Wi, CorPi;
    DMatrix Xi, Zscai, Zcori;
    IVector LinkWavei;
    
    getDatI(Y, Offset, Doffset, W, CorP, X, Zsca, Zcor, LinkWave, 
            I, J, Scuri, cor,
            Yi, Offseti, Doffseti, Wi, CorPi, Xi, Zscai, Zcori, LinkWavei);
    
    DVector PRi(s1), Vi(s1), V_Mui(s1); DMatrix Di(s1,pb);
    gee_prep(Yi, Xi, Offseti, LinkWavei, par, geestr, PRi, Di, Vi, V_Mui);
    DVector Phii(s1); DMatrix D2i(s1, pg);
    PhiandD2(LinkWavei, Doffseti, Zscai, par, geestr, Phii, D2i);
    DMatrix R(s1, s1), E(s2, pa);
    RandE(Zcori, CorPi, par, geestr, cor, R, E);
    
    HiandGi(PRi, Phii, Di, R, Vi, V_Mui, D2i, E, Wi, level, Hi, Gi);
    His(i) = Hi; Gis(i) = Gi;
  }
}

IVector genVI(IVector &Si, int c) {
  int s = Si.size(), k = 1;
  IVector ans(s * c); ans = 0;
  for (int i = 1; i <= s; i++) {
    if (Si(i) == 1) {
      for (int j = 1; j <= c; j++) {
        ans(k) = 1;
        k++;
      }
    }
  }
  return ans;
}

IVector genCrossVI(IVector &Si, int c) {
  int s = Si.size();
  IVector ans(s * (s - 1) * c * c / 2); ans = 0;
  IVector vv(c * c); vv = 1;
  Index1D I(0,0);
  for (int i = 1; i <= s - 1; i++) {
    for (int j = i + 1; j <= s; j++) {
      I = Index1D(1, c * c) + I.ubound();
      if (Si(i) == 1 && Si(j) == 1) 
        VecSubs(ans, I) = vv;
    }
  }
  return ans;
}

/*********************************************************/

DMatrix asDMatrix(SEXP a) {
  double *x;
  x = NUMERIC_POINTER(AS_NUMERIC(a));
  int *dims = INTEGER_POINTER(AS_INTEGER(GET_DIM(a)));
  DMatrix ans(dims[0], dims[1], x);
  return ans;
}

DVector asDVector(SEXP a) {
  double *x;
  x = NUMERIC_POINTER(AS_NUMERIC(a));
  int len = GET_LENGTH(a);
  DVector ans(len, x);
  return ans;
}

IVector asIVector(SEXP a) {
  int *x;
  x = INTEGER_POINTER(AS_INTEGER(a));
  int len = GET_LENGTH(a);
  IVector ans(len, x);
  return ans;
}

Vector<DVector> asVDVector(SEXP a) {//a is a matrix
  double *x;
  x = NUMERIC_POINTER(AS_NUMERIC(a));
  int *dims = INTEGER_POINTER(AS_INTEGER(GET_DIM(a)));
  Vector<DVector> ans(dims[1]);
  for (int i = 1; i <= ans.size(); i++) {
    DVector tmp(dims[0], x);
    ans(i) = tmp;
    x += dims[0];
  }
  return ans;
}

SEXP asSEXP(const DMatrix &a) {
  int size = a.num_cols() * a.num_rows();
  
  SEXP val;
  PROTECT(val = NEW_NUMERIC(size));
  double *p = NUMERIC_POINTER(val);
  const double *q = a.begin();
  for (int i = 0; i < size; i++) p[i] = q[i];
  //  SET_CLASS(val, ScalarString(mkChar("matrix")));
  
  SEXP dim;
  PROTECT(dim = NEW_INTEGER(2));
  INTEGER(dim)[0] = a.num_rows(); INTEGER(dim)[1] = a.num_cols();
  SET_DIM(val, dim);
  
  UNPROTECT(2);
  return val;
}

SEXP asSEXP(const DVector &a) {
  int size = a.size();
  SEXP val;
  PROTECT(val = NEW_NUMERIC(size));
  double *p = NUMERIC_POINTER(val);
  const double *q = a.begin();
  for (int i = 0; i < size; i++) p[i] = q[i];
  //  SET_CLASS(val, ScalarString(mkChar("vector")));
  
  SEXP len;
  PROTECT(len = NEW_INTEGER(1));
  INTEGER(len)[0] = size;
  SET_LENGTH(val, size);
  UNPROTECT(2);
  return val;
}

SEXP asSEXP(const IVector &a) {
  int size = a.size();
  SEXP val;
  PROTECT(val = NEW_INTEGER(size));
  int *p = INTEGER_POINTER(val);
  const int *q = a.begin();
  for (int i = 0; i < size; i++) p[i] = q[i];
  //  SET_CLASS(val, ScalarString(mkChar("vector")));
  
  SEXP len;
  PROTECT(len = NEW_INTEGER(1));
  INTEGER(len)[0] = size;
  SET_LENGTH(val, size);
  UNPROTECT(2);
  return val;
}


Control asControl(SEXP con) {
  //con is a list of trace, jack, j1s, fij, maxiter, epsilon
  int trace, jack, j1s, fij, maxiter;
  double tol;
  trace = INTEGER(VECTOR_ELT(con, 0))[0];
  jack = INTEGER(VECTOR_ELT(con, 1))[0];
  j1s = INTEGER(VECTOR_ELT(con, 2))[0];
  fij = INTEGER(VECTOR_ELT(con, 3))[0];
  maxiter = INTEGER(VECTOR_ELT(con, 4))[0];
  tol = REAL(VECTOR_ELT(con, 5))[0];
  Control Con(trace, jack, j1s, fij, maxiter, tol);
  return Con;
}

GeeParam asGeeParam(SEXP par) {
  //par is a list of beta, alpha, gamma;
  DVector Beta = asDVector(VECTOR_ELT(par, 0));
  DVector Alpha = asDVector(VECTOR_ELT(par, 1));
  DVector Gamma = asDVector(VECTOR_ELT(par, 2));
  GeeParam Par(Beta, Alpha, Gamma);
  return Par;
}

GeeStr asGeeStr(SEXP geestr) {
  //geestr is a list of maxwave, meanlink, v, scalelink, corrlink, scale.fix;
  int maxwave = INTEGER(AS_INTEGER(VECTOR_ELT(geestr, 0)))[0];
  IVector MeanLink = asIVector(AS_INTEGER(VECTOR_ELT(geestr, 1)));
  IVector V = asIVector(AS_INTEGER(VECTOR_ELT(geestr, 2)));
  IVector ScaleLink = asIVector(AS_INTEGER(VECTOR_ELT(geestr, 3)));
  int corrlink = INTEGER(AS_INTEGER(VECTOR_ELT(geestr, 4)))[0];
  int scalefix = INTEGER(AS_INTEGER(VECTOR_ELT(geestr, 5)))[0];
  GeeStr G(maxwave, MeanLink, V, ScaleLink, corrlink, scalefix);
  return G;
}

Corr asCorr(SEXP cor) {
  //cor is a list of corst, maxwave
  int corstr, maxwave;
  corstr = INTEGER(VECTOR_ELT(cor, 0))[0];
  maxwave = INTEGER(VECTOR_ELT(cor, 1))[0];
  Corr Cor(corstr, maxwave);
  return Cor;
}

SEXP asSEXP(GeeParam &Par) {
  SEXP ans;
  PROTECT(ans = NEW_LIST(19));
  SET_VECTOR_ELT(ans, 0, asSEXP(Par.beta()));
  SET_VECTOR_ELT(ans, 1, asSEXP(Par.alpha()));
  SET_VECTOR_ELT(ans, 2, asSEXP(Par.gamma()));
  SET_VECTOR_ELT(ans, 3, asSEXP(Par.vbeta()));
  SET_VECTOR_ELT(ans, 4, asSEXP(Par.valpha()));
  SET_VECTOR_ELT(ans, 5, asSEXP(Par.vgamma()));
  SET_VECTOR_ELT(ans, 6, asSEXP(Par.vbeta_naiv()));
  SET_VECTOR_ELT(ans, 7, asSEXP(Par.valpha_naiv()));
  SET_VECTOR_ELT(ans, 8, asSEXP(Par.valpha_stab()));
  SET_VECTOR_ELT(ans, 9, asSEXP(Par.vbeta_ajs()));
  SET_VECTOR_ELT(ans, 10, asSEXP(Par.valpha_ajs()));
  SET_VECTOR_ELT(ans, 11, asSEXP(Par.vgamma_ajs()));
  SET_VECTOR_ELT(ans, 12, asSEXP(Par.vbeta_j1s()));
  SET_VECTOR_ELT(ans, 13, asSEXP(Par.valpha_j1s()));
  SET_VECTOR_ELT(ans, 14, asSEXP(Par.vgamma_j1s()));
  SET_VECTOR_ELT(ans, 15, asSEXP(Par.vbeta_fij()));
  SET_VECTOR_ELT(ans, 16, asSEXP(Par.valpha_fij()));
  SET_VECTOR_ELT(ans, 17, asSEXP(Par.vgamma_fij()));
  
  IVector Err(1); Err(1) = Par.err();
  SET_VECTOR_ELT(ans, 18, asSEXP(Err));
  UNPROTECT(1);
  return ans;
}

/*********************************************************/

//class Control
Control::Control(int trace, int ajs, int j1s, int fij, 
                 int maxiter, double tol) :
  _trace(trace), _ajs(ajs), _j1s(j1s), _fij(fij),
  _maxiter(maxiter), _tol(tol){}
Control::Control(int *con, double tol) {
  _trace = con[0]; _ajs = con[1]; _j1s = con[2]; _fij = con[3];
  _maxiter = con[4]; _tol = tol;
}
Control::Control(const Control &C) :
  //{
  _trace(C.trace()), _ajs(C.ajs()), _j1s(C.j1s()),
  _fij(C.fij()), _maxiter(C.maxiter()), _tol(C.tol()) {}
// _trace = C.trace(); _ajs = C.ajs(); _j1s = C.j1s();
//_fij = C.fij(); _maxiter = C.maxiter(); _tol = C.tol();
//}

//class GeeParam
GeeParam::GeeParam(DVector Beta, DVector Alpha, DVector Gamma):
  _beta(Beta), _alpha(Alpha), _gamma(Gamma), _err(0) {
  int p = Beta.size(), q = Alpha.size(), r = Gamma.size();
  DMatrix vb(p,p), va(q,q), vg(r,r);
  _vbeta = vb; _vbeta_naiv = vb; _vbeta_ajs = vb; _vbeta_j1s = vb; _vbeta_fij = vb;
  _valpha = va; _valpha_naiv = va; _valpha_ajs = va; _valpha_j1s = va; _valpha_fij = va; _valpha_stab = va;
  _vgamma = vg; _vgamma_ajs = vg; _vgamma_j1s = vg; _vgamma_fij = vg;
}

GeeParam::GeeParam(DVector Beta, DVector Alpha, DVector Gamma,
                   DMatrix VBeta, DMatrix VBeta_naiv, 
                   DMatrix VBeta_ajs, DMatrix VBeta_j1s, 
                   DMatrix VBeta_fij,
                   DMatrix VAlpha, DMatrix VAlpha_stab,
                   DMatrix VAlpha_naiv, DMatrix VAlpha_ajs, 
                   DMatrix VAlpha_j1s, DMatrix VAlpha_fij,
                   DMatrix VGamma, DMatrix VGamma_ajs, 
                   DMatrix VGamma_j1s, DMatrix VGamma_fij):
  _beta(Beta),
  _alpha(Alpha),
  _gamma(Gamma),
  _vbeta(VBeta), _vbeta_naiv(VBeta_naiv),
  _vbeta_ajs(VBeta_ajs), _vbeta_j1s(VBeta_j1s),
  _vbeta_fij(VBeta_fij),
  _valpha(VAlpha), _valpha_stab(VAlpha_stab), 
  _valpha_naiv(VAlpha_naiv), _valpha_ajs(VAlpha_ajs), 
  _valpha_j1s(VAlpha_j1s), _valpha_fij(VAlpha_fij),
  _vgamma(VGamma),
  _vgamma_ajs(VGamma_ajs), _vgamma_j1s(VGamma_j1s),
  _vgamma_fij(VGamma_fij) {}

/*********************************************************/

void VecPrint(const DVector &v) {
  for (int i = 0; i < v.dim(); i++) Rprintf("%f ", v[i]);
  Rprintf("\n");
}

Fortran_Matrix<double> ident (int n) {
  Fortran_Matrix<double> ans(n,n);
  for (int i = 1; i <= n; i++) ans(i,i) = 1.0;
  return ans;
}

Fortran_Matrix<double> MatRowCol(const Fortran_Matrix<double> &mat, const Vector<double> &r, const Vector<double> &c) {
  int m = r.size(), n = c.size();
  Fortran_Matrix<double> ans(m,n);
  for (int i = 1; i <= m; i++)
    for (int j = 1; j <= n; j++) 
      ans(i,j) = mat((int) r(i), (int) c(j));
  return ans;
}
Fortran_Matrix<double> rho2mat(const Vector<double> &rho) {
  int s = rho.size(); // s = n(n-1)/2
  int n = (int) (0.5 * ( 1 + sqrt(1.0 + 8 * s)));
  Fortran_Matrix<double> fullmat = ident(n); 
  int k = 1;
  for (int i = 1; i <= n - 1; i++)
    for (int j = i + 1; j <= n; j++) {
      fullmat(i, j) = rho(k++);
      fullmat(j, i) = fullmat(i, j);
    }
    return fullmat;
}

DMatrix solve(const DMatrix &a, const DMatrix &b) {
  Subscript m = a.dim(1); // assert(m == a.dim(2));
  Subscript n = b.dim(1); // assert(m == n);
  Subscript l = b.dim(2);
  Vector<Subscript> index(m);
  DMatrix T(a), B(b);
  DMatrix ans(n,l);
  if (LU_factor(T, index) != 0) {
    // cerr << "LU_factor() failed." << endl; 
    return ans;
  }
  DVector v(m);
  for (int i  = 1; i <= l; i++) {
    v = asVec(MatCol(B,i));
    LU_solve(T, index, v);
    MatCol(ans, i) = asColMat(v);
  }
  return ans;  
}

DVector solve(const DMatrix &A, const DVector &b) {
  DMatrix T(A); Vector<Subscript> index(b.size());
  DVector ans(b);
  if (LU_factor(T, index) !=0) {
    //cerr << "LU_factor() failed." << endl;
    return ans;
  }
  
  if (LU_solve(T, index, ans) != 0)  {
    //cerr << "LU_Solve() failed." << endl;
    return ans;
  }
  return ans;
} 

DMatrix solve(const DMatrix &a) {
  DMatrix b = ident(a.dim(1));
  return solve(a, b);
}

DMatrix AtBiC(const DMatrix &A, const DMatrix &B, const DMatrix &C) {
  DMatrix BiC = solve(B, C);
  return Transpose_view(A) * BiC;
}

DVector AtBiC(const DMatrix &A, const DMatrix &B, const DVector &C) {
  DVector BiC = solve(B, C);
  return Transpose_view(A) * BiC;
}


DMatrix apply_elwise(const DMatrix &x, double f(double)) {
  DMatrix ans = x;
  for (int i = 1; i <= x.dim(1); i++)
    for (int j = 1; j <= x.dim(2); j++)
      ans(i, j) = f(x(i, j));
  return ans;
}

DVector apply_elwise(const DVector &x, double f(double)) {
  DVector ans = x;
  for (int i = 1; i <= x.dim(); i++)
    ans(i) = f(x(i));
  return ans;
}

DVector sqrt(const DVector &x) {
  return apply_elwise(x, sqrt);
}

double square(double x) {
  return x * x;
}

DVector square(const DVector &x) {
  return apply_elwise(x, square);
}

double reciproot(double x) {
  return 1./sqrt(x);
}


DVector reciproot(const DVector &x) {
  return apply_elwise(x, reciproot);
}

double recip(double x) {return 1./x;}

DVector recip(const DVector &x) {
  return apply_elwise(x, recip);
}

int cluscount(DVector &ID) {
  int ans = 1;
  for (int i = 1; i < ID.dim(); i++)
    if (ID(i - 1) != ID(i)) ans++;
    return ans;
}

Vector<int> clussize(DVector &ID) {
  int K = ID.size();
  Vector<int> ans(K); 
  ans = 1;
  //double id = ID(0);
  int k = 1;
  for (int i = 1; i <= (ID.dim() - 1); i++) {
    if (ID(i + 1) == ID(i)) ans(k) += 1;
    else k++;
  }
  return ans;
}

DVector SMult(const DVector &v1, const DVector &v2) {
  // assert (v1.dim() == v2.dim());
  DVector ans = v1;
  for (int i = 1; i <= v1.dim(); i++)
    ans(i) = v1(i) * v2(i);
  return ans;
}

DMatrix SMult(const DVector &v, const DMatrix &m) {
  // assert (v.dim() == m.dim(1));
  DMatrix tmp = m;
  for (int i = 1; i <= m.dim(1); i++)
    for (int j = 1; j <= m.dim(2); j++)
      tmp(i, j) = v(i) * m(i, j);
  return tmp;
}

DMatrix operator*(const DVector &v, const DMatrix &m) {
  return SMult(v, m);
}

DMatrix diag(const DVector &v) {
  int n = v.dim();
  DMatrix ans(n, n); ans = .0;
  for (int i = 1; i <= n; i++) ans(i, i) = v(i);
  return ans;
}

DVector diag(const DMatrix &m) {
  int n = m.dim(1); //assert m.dim(0) == m.dim(1);
  DVector ans(n); ans = .0;
  for (int i = 1; i <= n; i++) ans(i) = m(i,i);
  return ans;
}

DMatrix inv(const DMatrix &x) {
  return solve(x);
}

DMatrix fabs(const DMatrix &m) {
  return apply_elwise(m, fabs);
}

DVector fabs(const DVector &v) {
  return apply_elwise(v, fabs);
}

/*********************************************************/

// [[Rcpp::export]]
Rcpp::List getHessGrads(SEXP y, SEXP x, SEXP offset, SEXP doffset, SEXP w,
                       SEXP linkwave, SEXP zsca, SEXP zcor, SEXP corp,
                       SEXP clusz, SEXP geestr, SEXP cor, SEXP par) {
  DVector Y = asDVector(y), Offset = asDVector(offset), Doffset = asDVector(doffset), W = asDVector(w);
  IVector LinkWave = asIVector(linkwave); 
  DVector CorP = asDVector(corp); 
  DMatrix X = asDMatrix(x), Zsca = asDMatrix(zsca), Zcor = asDMatrix(zcor);
  IVector Clusz = asIVector(clusz); 
  //Control Con = asControl(con);   
  GeeParam Par = asGeeParam(par);   
  GeeStr Geestr = asGeeStr(geestr);   
  Corr Cor = asCorr(cor);   
  IVector Scur(Y.size()); Scur = 1;
  
  IVector level(2); level = 0;
  if (Geestr.ScaleFix() != 1) level(1) = 1;
  if (Cor.nparam() > 0) level(2) = 1;
  
  int K = Clusz.size();
  IVector Jack(K), ZcorSize(K);
  if (Cor.corst() > AR1) // == UNSTRUCTRUED || USERDEFINED || FIXED
    for (int i = 1; i <= K; i++) 
      ZcorSize(i) = Clusz(i) * (Clusz(i) - 1) / 2;
  else ZcorSize = 1;
  
  Hess Hi(Par), H(Par); Grad Gi(Par);
  Rcpp::List ans(6 + 3*Clusz.size());

  int pb = Par.p(), pa = Par.q(), pg = Par.r();
  
  //output
  Index1D I(0,0), J(0,0);
  int N = Clusz.size();
  DVector V0(pb);
  //cout << "N = " << N;
  for (int i = 1; i <= N; i++) {
    int s1 = Clusz(i), s2 = ZcorSize(i);
    I = Index1D(1, s1) + I.ubound();
    if (s2 > 0) J = Index1D(1, s2) + J.ubound();
    
    IVector Scuri = asVec(VecSubs(Scur, I));
    if (sum(Scuri) == 0)  continue;
    //get and valid data i
    DVector Yi, Offseti, Doffseti, Wi, CorPi;
    DMatrix Xi, Zscai, Zcori;
    IVector LinkWavei;
    
    getDatI(Y, Offset, Doffset, W, CorP, X, Zsca, Zcor, LinkWave, 
            I, J, Scuri, Cor, Yi, Offseti, Doffseti, Wi, CorPi, Xi, Zscai, Zcori, LinkWavei);
    
    DVector PRi(s1), Vi(s1), V_Mui(s1); DMatrix Di(s1,pb);
    gee_prep(Yi, Xi, Offseti, LinkWavei, Par, Geestr, PRi, Di, Vi, V_Mui);
    DVector Phii(s1); DMatrix D2i(s1, pg);
    PhiandD2(LinkWavei, Doffseti, Zscai, Par, Geestr, Phii, D2i);
    DMatrix R(s1, s1), E(s2, pa);
    RandE(Zcori, CorPi, Par, Geestr, Cor, R, E);
    
    HiandGi(PRi, Phii, Di, R, Vi, V_Mui, D2i, E, Wi, level, Hi, Gi);
    H.inc(Hi); 
    ans(6 + 3*(i-1)) = Gi.U1();
    ans(6 + 3*(i-1)+1) = Gi.U2();
    ans(6 + 3*(i-1)+2) = Gi.U3();
  }
  
  ans(0) = H.A();
  ans(1) = H.B();
  ans(2) = H.C();
  ans(3) = H.D();
  ans(4) = H.E();
  ans(5) = H.F();
  return ans;
}
