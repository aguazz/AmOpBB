#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

// Declaration for BBK
arma::mat BBKCpp(double t, arma::vec x, arma::vec t_u, arma::mat bt_u,
                 double y, double T, arma::vec sigma);

// [[Rcpp::export]]
arma::mat BBKCpp(double t, arma::vec x, arma::vec t_u, arma::mat bt_u, 
                 double y, double T, arma::vec sigma) {
  
  // Number of sigmas and length of the time grid
  arma::uword n = sigma.n_elem;
  arma::uword N = t_u.n_elem;
  
  // Cast sigma, x and t_u to matrices
  arma::mat sigma_mat = arma::repmat(sigma, 1, N);
  arma::mat x_mat = arma::repmat(x, 1, N);
  arma::mat t_u_mat = arma::repmat(t_u.t(), n, 1);
  
  // Compute z_u
  arma::mat z = (x_mat % (T - t_u_mat) + y * (t_u_mat - t)) / (T - t);
  arma::mat nu = sigma_mat % sqrt((t_u_mat - t) / (T - t) % (T - t_u_mat));
  z = (bt_u - z) / nu;
  
  // Kernel
  arma::mat K = (x_mat - y) / (T - t) % (1 - arma::normcdf(z));
  K += nu / (T - t_u_mat) % arma::normpdf(z);
  
  return K;
  
}

// [[Rcpp::export]]
arma::mat bBBCpp(arma::vec sigma, arma::vec t = 0, double tol = 1e-3, 
                 double y = 0, double T = 1, int N = 2e2) {
  
  // How many boundaries to compute
  arma::uword n = sigma.n_elem;
  
  // Check if the discretization is given
  if (t.n_elem == 1) {
    
    t = arma::regspace(0, T / N, T); // Discretization
    
  } else {
    
    N = t.n_elem - 1; // Number of subintervals
    T = t(N); // Expiration date
    
  }
  
  // Differences (a vector of length N)
  arma::vec Delta = arma::diff(t);
  
  // Boundary
  arma::mat bnd(n, N + 1); 
  bnd.fill(y);
  bnd.col(N - 1) = sigma / 2 * sqrt(atan(1) * 4 * (T - t(N - 1)) / 2)  + y;
  
  // Boundary computation
  arma::mat K = arma::zeros(n, N + 1);
  double H;
  double aux;
  
  for (arma::uword i = N - 1; i > 0; i--) {
    
    arma::vec bnd_old = bnd.col(i);
    arma::vec t_u = t.subvec(i, N - 1); // Time from the (i + 1)-th element
    arma::mat bt_u = bnd.cols(i, N - 1); // Boundary from the (i + 1)-th element
    
    aux = 2 * (T - t(i - 1)) / ( (T - t(i - 1)) + (t(N - 1) - t(i - 1)) );
    // Evaluate the integral int_{t_{N - 1}}^{T} sqrt((u - t_i) / (T - u)) du
    H = ( T - t(i - 1) ) * 
      ( atan(1) * 4 / 2 - atan( sqrt( (t(N - 1) - t(i - 1)) / (T - t(N - 1)) ) ) ) +
      ( T - t(N - 1) ) * sqrt( (t(N - 1) - t(i - 1)) / (T - t(N - 1)) );
    
    // Fixed point algorithm
    double error = 1;
    while (error > tol) {
      
      // Evaluate the kernel
      K = BBKCpp(t(i - 1), bnd_old, t_u, bt_u, y, T, sigma);
      
      // Updating the boundary by approximating the integral
      // in the Volterra equation using the Riemman sum
      bnd.col(i - 1) = y + ( K * Delta.subvec(i - 1, N - 2) + 
                      sigma * H / ( 2 * sqrt(2 * atan(1) * 4 * (T - t(i -1))) ) ) * 
                      aux;
      
      // Relative error
      error = max(abs((bnd.col(i - 1) - bnd_old) / bnd.col(i - 1)));
      // Caution: the while will keep running until all errors associated to 
      // each of the sigma's are smaller than the tolerange
      
      // Update
      bnd_old = bnd.col(i - 1);
      
    }
    
  }
  
  return bnd;
  
}
