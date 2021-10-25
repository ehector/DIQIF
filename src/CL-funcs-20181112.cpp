// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;
using namespace RcppParallel;

struct logCLnormal_worker : public Worker
{
  // source
  const RVector<double> beta;
  const RMatrix<double> cov; 
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  
  // destination 
  double SUM;
  
  // initialize with source and destination
  logCLnormal_worker(const NumericVector& beta, const NumericMatrix& cov, const NumericMatrix& block_y, 
                     const NumericMatrix& block_x, const int& dim) 
    : beta(beta), cov(cov), block_y(block_y), block_x(block_x), dim(dim), SUM(0) {}
  
  logCLnormal_worker(const logCLnormal_worker& little_worker, Split) 
    : beta(little_worker.beta), cov(little_worker.cov), block_y(little_worker.block_y), block_x(little_worker.block_x), 
      dim(little_worker.dim), SUM(0) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    double sigma_2r = 0;
    double sigma_2t = 0;
    double rho_rt = 0;
    double mean_r = 0;
    double mean_t = 0;
    
    for(unsigned int i = begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=r+1; t < dim; t++){
          sigma_2r = cov(r,r);
          sigma_2t = cov(t,t);
          rho_rt = cov(r,t)/pow(sigma_2r*sigma_2t, 0.5);
          
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k = 0; k < beta.length(); ++k) {
            mean_r += as_scalar(block_x(it*dim+r, k)*beta[k]);
            mean_t += as_scalar(block_x(it*dim+t, k)*beta[k]);
          }
          
          SUM += (-log(pow(sigma_2r*sigma_2t*(1-pow(rho_rt,2.0)), 0.5)) - 
            (pow(block_y(r, it)-mean_r, 2.0)/sigma_2r + pow(block_y(t, it)-mean_t, 2.0)/sigma_2t - 
            2*rho_rt*(block_y(r, it)-mean_r)*(block_y(t, it)-mean_t)/pow(sigma_2r*sigma_2t,0.5))/
              (2*(1-pow(rho_rt, 2.0))));
        }
      }
    }
  }
  
  // join my SUM with that of another logCLnormal_worker
  void join(const logCLnormal_worker& rhs) { 
    SUM += rhs.SUM;
  }
};

struct logCLnormalCS_worker : public Worker
{
  // source
  const RVector<double> beta;
  const double sigma;
  const double rho;
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  
  // destination 
  double SUM;
  
  // initialize with source and destination
  logCLnormalCS_worker(const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, 
                     const NumericMatrix& block_x, const int& dim) 
    : beta(beta), sigma(sigma), rho(rho), block_y(block_y), block_x(block_x), dim(dim), SUM(0) {}
  
  logCLnormalCS_worker(const logCLnormalCS_worker& little_worker, Split) 
    : beta(little_worker.beta), sigma(little_worker.sigma), rho(little_worker.rho), block_y(little_worker.block_y), 
      block_x(little_worker.block_x), dim(little_worker.dim), SUM(0) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    double sigma_2 = pow(sigma, 2.0);
    double mean_r = 0;
    double mean_t = 0;
    
    for(unsigned int i = begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=r+1; t < dim; t++){
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k = 0; k < beta.length(); ++k) {
            mean_r += as_scalar(block_x(it*dim+r, k)*beta[k]);
            mean_t += as_scalar(block_x(it*dim+t, k)*beta[k]);
          }
          
          SUM += (-log(pow(sigma_2*sigma_2*(1-pow(rho,2.0)), 0.5)) - 
            (pow(block_y(r, it)-mean_r, 2.0)/sigma_2 + pow(block_y(t, it)-mean_t, 2.0)/sigma_2 - 
            2*rho*(block_y(r, it)-mean_r)*(block_y(t, it)-mean_t)/pow(sigma_2*sigma_2,0.5))/
              (2*(1-pow(rho, 2.0))));
        }
      }
    }
  }
  
  // join my SUM with that of another logCLnormal_worker
  void join(const logCLnormalCS_worker& rhs) { 
    SUM += rhs.SUM;
  }
};

struct logCLnormalAR1_worker : public Worker
{
  // source
  const RVector<double> beta;
  const double sigma;
  const double rho;
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  
  // destination 
  double SUM;
  
  // initialize with source and destination
  logCLnormalAR1_worker(const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, 
                     const NumericMatrix& block_x, const int& dim) 
    : beta(beta), sigma(sigma), rho(rho), block_y(block_y), block_x(block_x), dim(dim), SUM(0) {}
  
  logCLnormalAR1_worker(const logCLnormalAR1_worker& little_worker, Split) 
    : beta(little_worker.beta), sigma(little_worker.sigma), rho(little_worker.rho), block_y(little_worker.block_y), block_x(little_worker.block_x), 
      dim(little_worker.dim), SUM(0) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    double sigma_2 = pow(sigma, 2.0);
    double rho_rt = 0;
    double mean_r = 0;
    double mean_t = 0;
    
    for(unsigned int i = begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=r+1; t < dim; t++){
          rho_rt = pow(rho, t-r);
          
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k = 0; k < beta.length(); ++k) {
            mean_r += as_scalar(block_x(it*dim+r, k)*beta[k]);
            mean_t += as_scalar(block_x(it*dim+t, k)*beta[k]);
          }
          
          SUM += (-log(pow(sigma_2*sigma_2*(1-pow(rho_rt,2.0)), 0.5)) - 
            (pow(block_y(r, it)-mean_r, 2.0)/sigma_2 + pow(block_y(t, it)-mean_t, 2.0)/sigma_2 - 
            2*rho_rt*(block_y(r, it)-mean_r)*(block_y(t, it)-mean_t)/sigma_2)/
              (2*(1-pow(rho_rt, 2.0))));
        }
      }
    }
  }
  
  // join my SUM with that of another logCLnormal_worker
  void join(const logCLnormalAR1_worker& rhs) { 
    SUM += rhs.SUM;
  }
};

struct logCLnormalind_worker : public Worker
{
  // source
  const RVector<double> beta;
  const double sigma;
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  
  // destination 
  double SUM;
  
  // initialize with source and destination
  logCLnormalind_worker(const NumericVector& beta, const double& sigma, const NumericMatrix& block_y, 
                        const NumericMatrix& block_x, const int& dim) 
    : beta(beta), sigma(sigma), block_y(block_y), block_x(block_x), dim(dim), SUM(0) {}
  
  logCLnormalind_worker(const logCLnormalind_worker& little_worker, Split) 
    : beta(little_worker.beta), sigma(little_worker.sigma), block_y(little_worker.block_y), block_x(little_worker.block_x), 
      dim(little_worker.dim), SUM(0) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    double sigma_2 = pow(sigma, 2.0);
    double mean_r = 0;
    double mean_t = 0;
    
    for(unsigned int i = begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=r+1; t < dim; t++){
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k = 0; k < beta.length(); ++k) {
            mean_r += as_scalar(block_x(it*dim+r, k)*beta[k]);
            mean_t += as_scalar(block_x(it*dim+t, k)*beta[k]);
          }
          
          SUM += (-log(pow(sigma_2*sigma_2, 0.5)) - (pow(block_y(r, it)-mean_r, 2.0)/sigma_2 + pow(block_y(t, it)-mean_t, 2.0)/sigma_2)/2);
        }
      }
    }
  }
  
  // join my SUM with that of another logCLnormal_worker
  void join(const logCLnormalind_worker& rhs) { 
    SUM += rhs.SUM;
  }
};

struct eenormalmean_worker : public Worker
{
  // source
  const RVector<double> beta;
  const RMatrix<double> cov; 
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> SUM;
  
  // initialize with source and destination
  eenormalmean_worker(const NumericVector& beta, const NumericMatrix& cov, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                      const int& dim, const int& len, NumericMatrix SUM) 
    : beta(beta), cov(cov), block_y(block_y), block_x(block_x), dim(dim), len(len), SUM(SUM) {}
    
    void operator()(std::size_t begin, std::size_t end) {
      double sigma_2r = 0;
      double sigma_2t = 0;
      double rho_rt = 0;
      double mean_r = 0;
      double mean_t = 0;
      double SUM_1 = 0;
      
      for(unsigned int i=begin; i < end; i++) {
        int it = static_cast<int>(i);
        for(int r=0; r < (dim-1); r++) {
          for(int t=(r+1); t < dim; t++){
            sigma_2r = cov(r,r);
            sigma_2t = cov(t,t);
            rho_rt = cov(r,t)/pow(sigma_2r*sigma_2t, 0.5);
            
            mean_r = 0;
            mean_t = 0;
            for (unsigned int k=0; k < beta.size(); k++){
              mean_r += block_x(it*dim+r, k)*beta[k];
              mean_t += block_x(it*dim+t, k)*beta[k];
            }
            
            for(unsigned int k=0; k < beta.size(); k++){
              SUM_1 = 0;
              for (unsigned int q=0; q < beta.size(); q++){
                SUM_1 += (block_x(it*dim+r, k)*block_x(it*dim+t, q) + block_x(it*dim+t, k)*block_x(it*dim+r, q))*beta[q];
              }
              SUM(k,it) = SUM(k,it) + 
                (block_x(it*dim+r,k)*(block_y(r, it)-mean_r)/(sigma_2r*(1-pow(rho_rt,2.0))) + 
                block_x(it*dim+t, k)*(block_y(t, it)-mean_t)/(sigma_2t*(1-pow(rho_rt,2.0))) - 
                rho_rt*(block_x(it*dim+t, k)*block_y(r, it)+block_x(it*dim+r,k)*block_y(t, it)-SUM_1)/
                  ((1-pow(rho_rt,2.0))*pow(sigma_2r*sigma_2t,0.5))); 
            }
          }
        }
      }
    }
};

struct eenormalderivmean_worker : public Worker
{
  // source
  const int p;
  const RMatrix<double> cov; 
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> SUM;
  
  // initialize with source and destination
  eenormalderivmean_worker(const int& p, const NumericMatrix& cov, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                      const int& dim, const int& len, NumericMatrix SUM) 
    : p(p), cov(cov), block_y(block_y), block_x(block_x), dim(dim), len(len), SUM(SUM) {}
  
  void operator()(std::size_t begin, std::size_t end) {

    double sigma_2r = 0;
    double sigma_2t = 0;
    double rho_rt = 0;
    int j = 0;
    
    for(unsigned int i=begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=(r+1); t < dim; t++){
          sigma_2r = cov(r,r);
          sigma_2t = cov(t,t);
          rho_rt = cov(r,t)/pow(sigma_2r*sigma_2t, 0.5);
          
          for (int k=0; k < p; k++){
            for (int q=0; q < p; q++){
              j = it*p + q;
              SUM(k,j) = SUM(k,j) + (sigma_2t*block_x(it*dim+r,k)*block_x(it*dim+r,q) + sigma_2r*block_x(it*dim+t,k)*block_x(it*dim+t,q)-
                rho_rt*pow(sigma_2r*sigma_2t, 0.5)*(block_x(it*dim+r,k)*block_x(it*dim+t,q)+block_x(it*dim+t,k)*block_x(it*dim+r,q)))/
                (sigma_2r*sigma_2t*(1-pow(rho_rt,2.0)));
            }
          }
        }
      }
    }
  }
};

struct psi_g_AR1derivvar_worker : public Worker
{
  // source
  const int p;
  const int d;
  const RVector<double> beta;
  const double sigma;
  const double rho;
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> SUM;
  
  // initialize with source and destination
  psi_g_AR1derivvar_worker(const int& p, const int&d, const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                           const int& dim, const int& len, NumericMatrix SUM) 
    : p(p), d(d), beta(beta), sigma(sigma), rho(rho), block_y(block_y), block_x(block_x), dim(dim), len(len), SUM(SUM) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    double sigma_2 = pow(sigma, 2.0);
    double sigma_3 = pow(sigma, 3.0);
    double rho_2tr;
    double rho_2r;
    double rho_2t;
    double rho_t_r;
    double mean_r;
    double mean_t;
    double cent_y_r;
    double cent_y_t;
    double z;
    double x_r_k;
    double x_t_k;
    double y_t_i;
    double y_r_i;
    double SUM_1 = 0;
    double SUM_2 = 0;
    int j = 0;
    
    for(unsigned int i=begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=(r+1); t < dim; t++){
          
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k=0; k < beta.size(); k++){
            mean_r += block_x(it*dim+r, k)*beta[k];
            mean_t += block_x(it*dim+t, k)*beta[k];
          }
          rho_2tr = pow(rho, 2.0*(t-r));
          rho_2r = pow(rho, 2*r);
          rho_2t = pow(rho, 2*t);
          rho_t_r = pow(rho,t-r);
          
          y_t_i = block_y(t, i);
          y_r_i = block_y(r, i);
          cent_y_r = y_r_i - mean_r;
          cent_y_t = y_t_i - mean_t;
          
          z = pow(cent_y_r, 2.0) + pow(cent_y_t, 2.0) - 2*rho_t_r*(cent_y_r)*(cent_y_t);
          
          for (int k=0; k < p; k++){
            x_r_k = block_x(it*dim+r, k);
            x_t_k = block_x(it*dim+t, k);
            
            SUM_1 = 0;
            for (unsigned int q=0; q < beta.size(); q++){
              SUM_1 += (x_r_k*block_x(it*dim+t, q) + x_t_k*block_x(it*dim+r, q))*beta[q];
            }
            SUM_2 = x_r_k*y_t_i + x_t_k*y_r_i - SUM_1;
            
            /* Update derivatives of psi_jk:*/
            // with respect to beta
            
            for (int q=0; q < p; q++){
              j = it*(p+d) + q;
              SUM(k,j) = SUM(k,j) +( (-x_r_k*block_x(it*dim+r,q) - x_t_k*block_x(it*dim+t,q) +
                rho_t_r*( x_r_k*block_x(it*dim+t,q) + x_t_k*block_x(it*dim+r,q) ))/
                  ( sigma_2*( 1-rho_2tr ) ));
            }
            // with respect to sigma
            SUM(k, it*(p+d)+p) = SUM(k, it*(p+d)+p) + (2* 
              ( rho_t_r* SUM_2 - x_r_k*cent_y_r - x_t_k*cent_y_t ) / (sigma_3*( 1-rho_2tr )));
            // with respect to rho
            SUM(k, it*(p+d)+p+1) = SUM(k, it*(p+d)+p+1) + ( 
              (r-t)*pow(rho, r+t-1)*(SUM_2*(pow(rho,2*r)+pow(rho,2*t))-2*(x_r_k*cent_y_r+x_t_k*cent_y_t)*pow(rho,r+t))/
                (sigma_2*pow(rho_2r-rho_2t,2.0))
              );
            
            /* Update derivatives of g_jk(sigma):*/
            // with respect to beta
            SUM(p, it*(p+d)+k) = SUM(p, it*(p+d)+k) + 2*( 
              -x_r_k*(cent_y_r) - x_t_k*(cent_y_t)
              + rho_t_r*SUM_2
            )/(sigma_3*(1-rho_2tr));
            
            /* Update derivates of g_jk(rho):*/
            // with respect to beta
            SUM(p+1, it*(p+d)+k) = SUM(p+1, it*(p+d)+k) + (t-r)*pow(rho, 2*(t-r)-1)/(pow(1-rho_2tr, 2.0)*sigma_2)*
            (2*(x_r_k*cent_y_r + x_t_k*cent_y_t - pow(rho, t-r)*SUM_2))
            - SUM_2*(t-r)*pow(rho,t-r-1)/(sigma_2*(1-rho_2tr));
            
          }
          
          /* Update derivates of g_jk(sigma):*/
          // with respect to sigma
          SUM(p, it*(p+d)+p) = SUM(p, it*(p+d)+p) +( 2/sigma_2 - 3*z/(pow(sigma, 4.0)*(1-rho_2tr)));
          
          // with respect to rho
          SUM(p, it*(p+d)+p+1) = SUM(p, it*(p+d)+p+1) + 2*(r-t)*
            (
                cent_y_r*cent_y_t * (pow(rho, t-r-1) + pow(rho, 3*(t-r)-1)) - 
               (pow(cent_y_r,2.0)+pow(cent_y_t,2.0))*pow(rho, 2*(t-r)-1)
            )/
              (sigma_3*pow(1-rho_2tr,2.0));
          
          /* Update derivates of g_jk(rho):*/
          // with respect to sigma
          SUM(p+1, it*(p+d)+p) = SUM(p+1, it*(p+d)+p) + 2*(t-r)*pow(rho,(t-r)-1)/(sigma_3*(1-rho_2tr))*(
            z*pow(rho, t-r)/(1-rho_2tr) - (cent_y_r)*(cent_y_t)
          )
          ;
          
          // with respect to rho
          if((t-r)<150){
            SUM(p+1, it*(p+d)+p+1) = SUM(p+1, it*(p+d)+p+1) + (
              (t-r)*(1-(2*r-2*t+1)*pow(rho,2*(r-t)))/pow(rho-pow(rho,2*(r-t)+1),2.0) +
                (z/sigma_2)*(r-t) * ((2*r-2*t+1)*pow(rho, 4*(r-t)-2) + (2*r-2*t-1)*pow(rho, 2*(r-t-1))) / pow(1-pow(rho, 2*(r-t)),3.0) +
                2*pow(r-t,2.0) * pow(rho,3*(t-r)-2) * cent_y_r*cent_y_t / (sigma_2*pow(1-rho_2tr, 2.0)) +
                (r-t)*cent_y_r*cent_y_t * ((r-t+1)*pow(rho, t-r-2)+(r-t-1)*pow(rho, 3*(t-r)-2))/
                  (sigma_2*pow(1-rho_2tr,2.0)) ); 
          }
          
        }
      }
    }
  }
};

struct psi_g_CSderivvar_worker : public Worker
{
  // source
  const int p;
  const int d;
  const RVector<double> beta;
  const double sigma;
  const double rho;
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> SUM;
  
  // initialize with source and destination
  psi_g_CSderivvar_worker(const int& p, const int&d, const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                           const int& dim, const int& len, NumericMatrix SUM) 
    : p(p), d(d), beta(beta), sigma(sigma), rho(rho), block_y(block_y), block_x(block_x), dim(dim), len(len), SUM(SUM) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    double sigma_2 = pow(sigma, 2.0);
    double sigma_3 = pow(sigma, 3.0);
    double rho_2 = pow(rho, 2.0);
    double mean_r;
    double mean_t;
    double cent_y_r;
    double cent_y_t;
    double z;
    double x_r_k;
    double x_t_k;
    double y_t_i;
    double y_r_i;
    double SUM_1 = 0;
    double SUM_2 = 0;
    int j = 0;
    
    for(unsigned int i=begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=(r+1); t < dim; t++){
          
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k=0; k < beta.size(); k++){
            mean_r += block_x(it*dim+r, k)*beta[k];
            mean_t += block_x(it*dim+t, k)*beta[k];
          }
          
          y_t_i = block_y(t, i);
          y_r_i = block_y(r, i);
          cent_y_r = y_r_i - mean_r;
          cent_y_t = y_t_i - mean_t;
          
          z = pow(cent_y_r, 2.0) + pow(cent_y_t, 2.0) - 2*rho*(cent_y_r)*(cent_y_t);
          
          for (int k=0; k < p; k++){
            x_r_k = block_x(it*dim+r, k);
            x_t_k = block_x(it*dim+t, k);
            
            SUM_1 = 0;
            for (unsigned int q=0; q < beta.size(); q++){
              SUM_1 += (x_r_k*block_x(it*dim+t, q) + x_t_k*block_x(it*dim+r, q))*beta[q];
            }
            SUM_2 = x_r_k*y_t_i + x_t_k*y_r_i - SUM_1;
            
            /* Update derivatives of psi_jk:*/
            // with respect to beta
            for (int q=0; q < p; q++){
              j = it*(p+d) + q;
              SUM(k,j) = SUM(k,j) +( (-x_r_k*block_x(it*dim+r,q) - x_t_k*block_x(it*dim+t,q) +
                rho*( x_r_k*block_x(it*dim+t,q) + x_t_k*block_x(it*dim+r,q) ))/( sigma_2*( 1-rho_2 ) )
              );
            }
            // with respect to sigma
            SUM(k, it*(p+d)+p) = SUM(k, it*(p+d)+p) + (2* 
              ( rho* SUM_2 - x_r_k*cent_y_r - x_t_k*cent_y_t ) / (sigma_3*( 1-rho_2 )));
            // with respect to rho
            SUM(k, it*(p+d)+p+1) = SUM(k, it*(p+d)+p+1) + 1/(sigma_2*pow(1-rho_2,2.0))*(
            2*rho*(x_r_k*cent_y_r+x_t_k*cent_y_t) - (1+rho_2)*SUM_2
            );
            
            /* Update derivatives of g_jk(sigma):*/
            // with respect to beta
            SUM(p, it*(p+d)+k) = SUM(p, it*(p+d)+k) + 2*( 
              -x_r_k*(cent_y_r) - x_t_k*(cent_y_t)
              + rho*SUM_2)/(sigma_3*(1-rho_2)
            );
            
            /* Update derivatives of g_jk(rho):*/
            // with respect to beta
            SUM(p+1, it*(p+d)+k) = SUM(p+1, it*(p+d)+k) + rho/(pow(1-rho_2, 2.0)*sigma_2)*
            (2*(x_r_k*cent_y_r + x_t_k*cent_y_t - rho*SUM_2))
              - SUM_2/(sigma_2*(1-rho_2)
            );
            
          }
          
          /* Update derivatives of g_jk(sigma):*/
          // with respect to sigma
          SUM(p, it*(p+d)+p) = SUM(p, it*(p+d)+p) +( 2/sigma_2 - 3*z/(pow(sigma, 4.0)*(1-rho_2)));
          
          // with respect to rho
          SUM(p, it*(p+d)+p+1) = SUM(p, it*(p+d)+p+1) + 2/(sigma_3*pow(1-rho_2,2.0))*(
            rho*(pow(cent_y_r,2.0) + pow(cent_y_t,2.0)) - (1+rho_2)*cent_y_r*cent_y_t
          );
          
          /* Update derivatives of g_jk(rho):*/
          // with respect to sigma
          SUM(p+1, it*(p+d)+p) = SUM(p+1, it*(p+d)+p) + 2/(sigma_3*(1-rho_2))*(
            z*rho/(1-rho_2) - (cent_y_r)*(cent_y_t)
          );
          
          // with respect to rho
          SUM(p+1, it*(p+d)+p+1) = SUM(p+1, it*(p+d)+p+1) + (
            (1+rho_2)/pow(1-rho_2,2.0) - (3*rho_2+1)*z/(sigma_2*pow(1-rho_2,3.0))+4*rho*(cent_y_r)*(cent_y_t)/(sigma_2*pow(1-rho_2,2.0))
            );
          
        }
      }
    }
  }
};

struct psi_g_indderivvar_worker : public Worker
{
  // source
  const int p;
  const int d;
  const RVector<double> beta;
  const double sigma;
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> SUM;
  
  // initialize with source and destination
  psi_g_indderivvar_worker(const int& p, const int&d, const NumericVector& beta, const double& sigma, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                          const int& dim, const int& len, NumericMatrix SUM) 
    : p(p), d(d), beta(beta), sigma(sigma), block_y(block_y), block_x(block_x), dim(dim), len(len), SUM(SUM) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    double sigma_2 = pow(sigma, 2.0);
    double sigma_3 = pow(sigma, 3.0);
    double mean_r;
    double mean_t;
    double cent_y_r;
    double cent_y_t;
    double z;
    double x_r_k;
    double x_t_k;
    double y_t_i;
    double y_r_i;
    double SUM_1 = 0;
    double SUM_2 = 0;
    int j = 0;
    
    for(unsigned int i=begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=(r+1); t < dim; t++){
          
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k=0; k < beta.size(); k++){
            mean_r += block_x(it*dim+r, k)*beta[k];
            mean_t += block_x(it*dim+t, k)*beta[k];
          }
          
          y_t_i = block_y(t, i);
          y_r_i = block_y(r, i);
          cent_y_r = y_r_i - mean_r;
          cent_y_t = y_t_i - mean_t;
          
          z = pow(cent_y_r, 2.0) + pow(cent_y_t, 2.0);
          
          for (int k=0; k < p; k++){
            x_r_k = block_x(it*dim+r, k);
            x_t_k = block_x(it*dim+t, k);
            
            SUM_1 = 0;
            for (unsigned int q=0; q < beta.size(); q++){
              SUM_1 += (x_r_k*block_x(it*dim+t, q) + x_t_k*block_x(it*dim+r, q))*beta[q];
            }
            SUM_2 = x_r_k*y_t_i + x_t_k*y_r_i - SUM_1;
            
            /* Update derivatives of psi_jk:*/
            // with respect to beta
            for (int q=0; q < p; q++){
              j = it*(p+d) + q;
              SUM(k,j) = SUM(k,j) +( (-x_r_k*block_x(it*dim+r,q) - x_t_k*block_x(it*dim+t,q))/( sigma_2 )
              );
            }
            // with respect to sigma
            SUM(k, it*(p+d)+p) = SUM(k, it*(p+d)+p) + (2* 
              ( - x_r_k*cent_y_r - x_t_k*cent_y_t ) / sigma_3);
            
            /* Update derivatives of g_jk(sigma):*/
            // with respect to beta
            SUM(p, it*(p+d)+k) = SUM(p, it*(p+d)+k) + 2*( 
              -x_r_k*(cent_y_r) - x_t_k*(cent_y_t) )/sigma_3;

          }
          
          /* Update derivatives of g_jk(sigma):*/
          // with respect to sigma
          SUM(p, it*(p+d)+p) = SUM(p, it*(p+d)+p) +( 2/sigma_2 - 3*z/(pow(sigma, 4.0)));
          
        }
      }
    }
  }
};


struct eenormalCSvar_worker : public Worker
{
  // source
  const RVector<double> beta;
  const double sigma;
  const double rho;
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> SUM;
  
  // initialize with source and destination
  eenormalCSvar_worker(const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                      const int& dim, const int& len, NumericMatrix SUM) 
    : beta(beta), sigma(sigma), rho(rho), block_y(block_y), block_x(block_x), dim(dim), len(len), SUM(SUM) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    double mean_r = 0;
    double mean_t = 0;
    double SUM_1 = 0;
    
    for(unsigned int i=begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=(r+1); t < dim; t++){
          
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k=0; k < beta.size(); k++){
            mean_r += block_x(it*dim+r, k)*beta[k];
            mean_t += block_x(it*dim+t, k)*beta[k];
          }
          
          for (unsigned int k=0; k < beta.size(); k++){
            SUM_1 = 0;
            for (unsigned int q=0; q < beta.size(); q++){
              SUM_1 += (block_x(it*dim+r, k)*block_x(it*dim+t, q) + block_x(it*dim+t, k)*block_x(it*dim+r, q))*beta[q];
            }
            
            SUM(k,it) = SUM(k, it) + (block_x(it*dim+r,k)*(block_y(r, i)-mean_r)/(pow(sigma,2.0)*(1-pow(rho,2.0))) + 
              block_x(it*dim+t,k)*(block_y(t, i)-mean_t)/(pow(sigma,2.0)*(1-pow(rho,2.0))) - 
              rho*(block_x(it*dim+t,k)*block_y(r, i)+block_x(it*dim+r,k)*block_y(t, i)-SUM_1)/((1-pow(rho,2.0))*pow(sigma,2.0)));
          }
          
          SUM(beta.size(),it) = SUM(beta.size(),it) + 
            as_scalar(-2/sigma + (pow(block_y(r, i)-mean_r, 2.0)-2*rho*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)+pow(block_y(t, i)-mean_t, 2.0))/
            ((1-pow(rho, 2.0))*pow(sigma,3.0)));
          
          SUM(beta.size()+1, it) = SUM(beta.size()+1, it) +
            as_scalar(rho/(1-pow(rho, 2.0))-(rho/pow(1-pow(rho,2.0),2.0))*(pow(block_y(r, i)-mean_r, 2.0)-2*rho*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)+
            pow(block_y(t, i)-mean_t, 2.0))/pow(sigma, 2.0)+((block_y(r, i)-mean_r)*(block_y(t, i)-mean_t))/(pow(sigma,2.0)*(1-pow(rho,2.0))));
        
        }
      }
    }
  }
};

struct eenormalAR1var_worker : public Worker
{
  // source
  const RVector<double> beta;
  const double sigma;
  const double rho;
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> SUM;
  
  // initialize with source and destination
  eenormalAR1var_worker(const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                       const int& dim, const int& len, NumericMatrix SUM) 
    : beta(beta), sigma(sigma), rho(rho), block_y(block_y), block_x(block_x), dim(dim), len(len), SUM(SUM) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    double mean_r = 0;
    double mean_t = 0;
    double SUM_1 = 0;
    
    for(unsigned int i=begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=(r+1); t < dim; t++){
          
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k=0; k < beta.size(); k++){
            mean_r += block_x(it*dim+r, k)*beta[k];
            mean_t += block_x(it*dim+t, k)*beta[k];
          }
          
          for (unsigned int k=0; k < beta.size(); k++){
            SUM_1 = 0;
            for (unsigned int q=0; q < beta.size(); q++){
              SUM_1 += (block_x(it*dim+r, k)*block_x(it*dim+t, q) + block_x(it*dim+t, k)*block_x(it*dim+r, q))*beta[q];
            }
            SUM(k,it) = SUM(k,it) +
              (block_x(it*dim+r,k)*(block_y(r, i)-mean_r)/(pow(sigma,2.0)*(1-pow(rho,2.0*(t-r)))) + 
              block_x(it*dim+t,k)*(block_y(t, i)-mean_t)/(pow(sigma,2.0)*(1-pow(rho,2.0*(t-r)))) - 
              pow(rho,t-r)*(block_x(it*dim+t,k)*block_y(r, i)+block_x(it*dim+r,k)*block_y(t, i)-SUM_1)/((1-pow(rho,2.0*(t-r)))*pow(sigma,2.0)));
          }
          
          SUM(beta.size(),it) = SUM(beta.size(),it) +
            as_scalar(-2/sigma + (pow(block_y(r, i)-mean_r, 2.0)-2*pow(rho,t-r)*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)+pow(block_y(t, i)-mean_t, 2.0))/
            ((1-pow(rho, 2.0*(t-r)))*pow(sigma,3.0)));
          
          SUM(beta.size()+1, it) = SUM(beta.size()+1, it) +
            as_scalar((r-t)/(rho-pow(rho, 2.0*(r-t)+1)) - (t-r)*pow(rho, 2.0*(t-r)-1)*
            (pow(block_y(r, i)-mean_r, 2.0)-2*pow(rho, t-r)*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)+pow(block_y(t, i)-mean_t,2.0))/
              (pow(1-pow(rho, 2.0*(t-r)), 2.0)*pow(sigma, 2.0)) +
                (t-r)*pow(rho, t-r-1)*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)/
                  (pow(sigma, 2.0)*(1-pow(rho, 2.0*(t-r)))));
        }
      }
    }
  }
};

struct eenormalindvar_worker : public Worker
{
  // source
  const RVector<double> beta;
  const double sigma;
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> SUM;
  
  // initialize with source and destination
  eenormalindvar_worker(const NumericVector& beta, const double& sigma, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                        const int& dim, const int& len, NumericMatrix SUM) 
    : beta(beta), sigma(sigma), block_y(block_y), block_x(block_x), dim(dim), len(len), SUM(SUM) {}
  
  void operator()(std::size_t begin, std::size_t end) {

    double mean_r = 0;
    double mean_t = 0;
    
    for(unsigned int i=begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=(r+1); t < dim; t++){
          
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k=0; k < beta.size(); k++){
            mean_r += block_x(it*dim+r, k)*beta[k];
            mean_t += block_x(it*dim+t, k)*beta[k];
          }
          
          for (unsigned int k=0; k < beta.size(); k++){
            SUM(k,it) = SUM(k,it) +
              (block_x(it*dim+r,k)*(block_y(r, i)-mean_r) + block_x(it*dim+t,k)*(block_y(t, i)-mean_t))/pow(sigma, 2.0);
          }
          
          SUM(beta.size(),it) = SUM(beta.size(),it) - 2/sigma + (pow(block_y(r, i)-mean_r, 2.0) + pow(block_y(t, i)-mean_t, 2.0))/pow(sigma, 3.0);
        }
      }
    }
  }
};

struct GEEnormalmean_worker : public Worker
{
  // source
  const RVector<double> beta;
  const RMatrix<double> cov_inv; 
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> SUM;
  
  // initialize with source and destination
  GEEnormalmean_worker(const NumericVector& beta, const NumericMatrix& cov_inv, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                      const int& dim, const int& len, NumericMatrix SUM) 
    : beta(beta), cov_inv(cov_inv), block_y(block_y), block_x(block_x), dim(dim), len(len), SUM(SUM) {}
  
  void operator()(std::size_t begin, std::size_t end) {
  
    for(unsigned int i=begin; i < end; i++) {
      double SUM_1 = 0;
      double SUM_2 = 0;
      double SUM_3 = 0;
      int it = static_cast<int>(i);
      
      for(unsigned int q=0; q<beta.length(); q++){
        SUM_3 = 0;
        for(unsigned int r=0; r<dim; r++){
          SUM_1 = 0;
          SUM_2 = 0;
          for(unsigned int s=0; s<beta.length(); s++){
            SUM_1 = SUM_1 + block_x(it*dim+r,s)*beta[s];
          }
          for(unsigned int s=0; s<dim; s++){
            SUM_2 = SUM_2 + block_x(it*dim+s,q)*cov_inv(s,r);
          }
          SUM_3 = SUM_3 + (block_y(r,it) - SUM_1)*SUM_2;
        }
        SUM(q,it) = SUM_3;
      }
    }
  }
};

struct GEEnormalderivmean_worker : public Worker
{
  // source
  const int p;
  const RMatrix<double> cov_inv; 
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> SUM;
  
  // initialize with source and destination
  GEEnormalderivmean_worker(const int& p, const NumericMatrix& cov_inv, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                            const int& dim, const int& len, NumericMatrix SUM) 
    : p(p), cov_inv(cov_inv), block_y(block_y), block_x(block_x), dim(dim), len(len), SUM(SUM) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    double SUM_1 = 0;
    int j = 0;
    
    for(unsigned int i=begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(unsigned int q=0; q<p; q++){
        for(unsigned int r=0; r<p; r++){
          j = it*p + r;
          SUM_1 = 0;
          for(unsigned int s=0; s<dim; s++){
            for(unsigned int t=0; t<dim; t++){
              SUM_1 = SUM_1 + block_x(it*dim+t,q)*cov_inv(t,s)*block_x(it*dim+s,r);
            }
          }
          SUM(q,j) = SUM(q,j) - SUM_1;
        }
      }
    }
  }
};



// [[Rcpp::export]]
double logCLnormal_par(const NumericVector& beta, const NumericMatrix& cov, const NumericMatrix& block_y, 
                       const NumericMatrix& block_x, const double& m, const double& n)
{
  int dim = m;
  int len = n;
  
  // logCLnormal functor (pass input matrices)
  logCLnormal_worker little_worker(beta, cov, block_y, block_x, dim);
  
  // call parallelReduce to do the work
  parallelReduce(0, len, little_worker);
  
  // return the SUM
  return little_worker.SUM;
}


// [[Rcpp::export]]
List eenormalmean_par(const NumericVector& beta, const NumericMatrix& cov, const NumericMatrix& block_y, 
                      const NumericMatrix& block_x, const double& m, const double& n)
{
  int dim = m;
  int len = n;
  NumericMatrix SUM (beta.length(), len);
  
  
  // eenormalmean_par functor (pass input and output matrices)
  eenormalmean_worker little_worker(beta, cov, block_y, block_x, dim, len, SUM);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  // return the output list
  List output(len);
  for (int i=0; i<len; i++){
    arma::mat MAT(SUM.column(i).begin(), beta.size(), 1, false);
    output[i] = MAT;
  }
  return output;
}

// [[Rcpp::export]]
List GEEnormalmean_par(const NumericVector& beta, const NumericMatrix& cov_inv, const NumericMatrix& block_y, 
                      const NumericMatrix& block_x, const double& m, const double& n)
{
  int dim = m;
  int len = n;
  NumericMatrix SUM (beta.length(), len);
  
  
  // eenormalmean_par functor (pass input and output matrices)
  GEEnormalmean_worker little_worker(beta, cov_inv, block_y, block_x, dim, len, SUM);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  // return the output list
  List output(len);
  for (int i=0; i<len; i++){
    arma::mat MAT(SUM.column(i).begin(), beta.size(), 1, false);
    output[i] = MAT;
  }
  return output;
}

// [[Rcpp::export]]
List GEEnormalderivmean_par(const NumericMatrix& cov_inv, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                           const double& m, const double& n)
{
  int dim = m;
  int len = n;
  int p = block_x.ncol();
  NumericMatrix SUM (p, p*len);
  
  
  // eenormalderivmean_par functor (pass input and output matrices)
  GEEnormalderivmean_worker little_worker(p, cov_inv, block_y, block_x, dim, len, SUM);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  // return the output list
  List output(len);
  for (int i=0; i<len; i++){
    arma::mat MAT(SUM.column(i*p).begin(), p, p, false);
    output[i] = MAT;
  }
  return output;
}

// [[Rcpp::export]]
List eenormalderivmean_par(const NumericMatrix& cov, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                           const double& m, const double& n)
{
  int dim = m;
  int len = n;
  int p = block_x.ncol();
  NumericMatrix SUM (p, p*len);
  
  
  // eenormalderivmean_par functor (pass input and output matrices)
  eenormalderivmean_worker little_worker(p, cov, block_y, block_x, dim, len, SUM);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  // return the output list
  List output(len);
  for (int i=0; i<len; i++){
      arma::mat MAT(SUM.column(i*p).begin(), p, p, false);
      output[i] = MAT;
  }
  return output;
}

// [[Rcpp::export]]
List psi_g_AR1derivvar_par(const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                           const double& m, const double& n)
{
  int dim = m;
  int len = n;
  int p = block_x.ncol();
  int d = 2;
  NumericMatrix SUM (p+d, (p+d)*len);
  
  
  // eenormalderivmean_par functor (pass input and output matrices)
  psi_g_AR1derivvar_worker little_worker(p, d, beta, sigma, rho, block_y, block_x, dim, len, SUM);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  List output(len);
  for (int i=0; i<len; i++){
    arma::mat MAT(SUM.column(i*(p+d)).begin(), p+d, p+d, false);
    output[i] = MAT;
  }
  return output;
}

// [[Rcpp::export]]
List psi_g_CSderivvar_par(const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                           const double& m, const double& n)
{
  int dim = m;
  int len = n;
  int p = block_x.ncol();
  int d = 2;
  NumericMatrix SUM (p+d, (p+d)*len);
  
  
  // eenormalderivmean_par functor (pass input and output matrices)
  psi_g_CSderivvar_worker little_worker(p, d, beta, sigma, rho, block_y, block_x, dim, len, SUM);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  List output(len);
  for (int i=0; i<len; i++){
    arma::mat MAT(SUM.column(i*(p+d)).begin(), p+d, p+d, false);
    output[i] = MAT;
  }
  return output;
}

// [[Rcpp::export]]
List psi_g_indderivvar_par(const NumericVector& beta, const double& sigma, const NumericMatrix& block_y, const NumericMatrix& block_x, 
                          const double& m, const double& n)
{
  int dim = m;
  int len = n;
  int p = block_x.ncol();
  int d = 1;
  NumericMatrix SUM (p+d, (p+d)*len);
  
  
  // eenormalderivmean_par functor (pass input and output matrices)
  psi_g_indderivvar_worker little_worker(p, d, beta, sigma, block_y, block_x, dim, len, SUM);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  List output(len);
  for (int i=0; i<len; i++){
    arma::mat MAT(SUM.column(i*(p+d)).begin(), p+d, p+d, false);
    output[i] = MAT;
  }
  return output;
}

// [[Rcpp::export]]
List eenormalCSvar_par(const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, 
                             const NumericMatrix& block_x, const double& m, const double& n)
{
  int dim = m;
  int len = n;
  int q = beta.size()+2;
  NumericMatrix SUM (q, len);
  
  
  // eenormalCSvar_par functor (pass input and output matrices)
  eenormalCSvar_worker little_worker(beta, sigma, rho, block_y, block_x, dim, len, SUM);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  // return the output list
  List output(len);
  for (int i=0; i<len; i++){
    arma::mat MAT(SUM.column(i).begin(), q, 1, false);
    output[i] = MAT;
  }
  return output;
}

// [[Rcpp::export]]
List eenormalAR1var_par(const NumericVector& beta, const double& sigma, const double& rho, const NumericMatrix& block_y, 
                        const NumericMatrix& block_x, const double& m, const double& n)
{
  int dim = m;
  int len = n;
  int q = beta.size()+2;
  NumericMatrix SUM (q, len);
  
  
  // eenormalAR1var_par functor (pass input and output matrices)
  eenormalAR1var_worker little_worker(beta, sigma, rho, block_y, block_x, dim, len, SUM);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  // return the output list
  List output(len);
  for (int i=0; i<len; i++){
    arma::mat MAT(SUM.column(i).begin(), q, 1, false);
    output[i] = MAT;
  }
  return output;
}

// [[Rcpp::export]]
List eenormalindvar_par(const NumericVector& beta, const double& sigma, const NumericMatrix& block_y, 
                        const NumericMatrix& block_x, const double& m, const double& n)
{
  int dim = m;
  int len = n;
  int q = beta.size()+1;
  NumericMatrix SUM (q, len);
  
  
  // eenormalindvar_par functor (pass input and output matrices)
  eenormalindvar_worker little_worker(beta, sigma, block_y, block_x, dim, len, SUM);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  // return the output list
  List output(len);
  for (int i=0; i<len; i++){
    arma::mat MAT(SUM.column(i).begin(), q, 1, false);
    output[i] = MAT;
  }
  return output;
}
