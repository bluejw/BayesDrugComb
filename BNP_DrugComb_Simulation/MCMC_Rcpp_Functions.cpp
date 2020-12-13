//[[Rcpp::depends(RcppArmadillo,RcppEigen)]]
# include <RcppArmadillo.h>
# include <RcppArmadilloExtensions/sample.h>
# include <RcppEigen.h> // warnings may be generated when include RcppEigen, but it will not affect the usage


using namespace std;
using namespace arma;
using namespace Rcpp; 


const double log2pi = log(2.0*M_PI);



// [[Rcpp::export]]
SEXP mat_mult_rcpp(Eigen::MatrixXd& A, Eigen::MatrixXd& B){
  
  // fast high-dimensional matrix multiplication
  //
  // args: A, B matrix
  // returns: C = A * B
  
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}



// [[Rcpp::export]]
double dmvn_rcpp(arma::rowvec& x, arma::rowvec& mean, arma::mat& sigma, bool logd = false){ 
  
  // calculate density of multivariate normal distribution
  //
  // args: x: row vector data
  //      mean: row vector mean, sigma: covaraicne matrix  
  //      logd: true for taking log
  // returns: out: pdf (or log pdf) of multivariate normal distribution
  
  int xdim = x.size(); 
  arma::mat rooti = trans(inv(trimatu(chol(sigma))));
  double rootisum = sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0)*log2pi;
  
  arma::vec z = rooti*trans(x-mean);
  double out = constants-0.5*sum(z%z)+rootisum;
  
  if (logd == false){ out = exp(out); }
  return(out);
}



// [[Rcpp::export]]
arma::mat rmvn_rcpp(const int n, arma::vec& mean, arma::mat& sigma){
  
  // randomly generate samples from multivariate normal distribution
  //
  // args: n: number of data 
  //      mean: row vector mean, sigma: covaraicne matrix  
  // returns: out: random samples from multivariate normal distribution
  
  int k = sigma.n_cols; // dimension of the multivariate normal distribution
  arma::mat z = randn(n, k);
  arma::mat out = repmat(mean,1,n).t()+z*chol(sigma);
  return(out);
}



// [[Rcpp::export]]
arma::mat riwish_rcpp(const int df, arma::mat& S){
  
  // randomly generate matrix from inverse wishart distribution
  //
  // args: df: degrees of freedom
  //       S: inverse of scale matrix 
  // returns: out: random matrix from inverse wishart distribution
  
  S = S.i(); 
  int m = S.n_rows;
  
  arma::mat Z(m,m); // Bartlett decomposition
  for (int i = 0; i < m; i++){
    Z(i,i) = sqrt(R::rchisq(df-i)); // Fill the diagonal
  }
  for (int j = 0; j < m; j++){  
    for(int i = j+1; i < m; i++){    
      Z(i,j) = R::rnorm(0,1); // Fill the lower matrix 
    }
  }
  
  arma::mat C = trimatl(Z).t() * chol(S);
  // Random matrix from inverse wishart distribution
  arma::mat out = (C.t()*C).i();
  return(out);
}



// [[Rcpp::export]]
double logpmf_pi_n_rcpp(arma::vec& pi_n, const double m_0, arma::uvec& sigma, arma::mat& Kappa){
  
  // log probability mass function of partition pi
  // 
  // args: pi_n: partition, m_0: mass parameter,  sigma: permutation,
  //       Kappa: ST kernel matrix, whose (i,j)-th element is the similarity between subject i and j 
  // returns: log_pmf: log probability mass function value of a partition pi_n  
  
  double log_pmf = 0; 
  double epsilon = 1e-300; 
  int n = pi_n.size(); // number of data in partition pi_n
  
  // pmf of partition pi is defined as 1 for t = 1
  if (n == 1) { return(log_pmf); }
  
  int current_subject; 
  double similarity_sum_same_class, similarity_sum;
  arma::vec Kappa_current_subject(n);
  
  // pmf of partition pi is defined iteratively for t > 1 
  for (int t=2; t<(n+1); t++){
    
    // current subject index at t-1: sigma_{t}
    current_subject = sigma(t-1)-1;
    // permutation at t-1: sigma_{1},...sigma_{t-1} 
    arma::uvec previous_subjects = sigma.subvec(0,t-2)-1;
    
    // previous partition at t-1: pi_n(sigma_{1},...sigma_{t-1})
    arma::vec previous_classes = pi_n.elem(previous_subjects); 
    // previous t-1 subjects in the same class of t-th subject 
    arma::uvec previous_subjects_same_class = previous_subjects.elem(find(previous_classes == pi_n[current_subject]));
    
    // if the class number of t-th subject 
    // exists in previous t-1 subjects' class numbers
    if (previous_subjects_same_class.size() > epsilon){
      
      // similarities between current subject and 
      // previous subjects in the same class or all previous subjects
      Kappa_current_subject = Kappa.col(current_subject);
      similarity_sum_same_class = sum(Kappa_current_subject.elem(previous_subjects_same_class));
      similarity_sum = sum(Kappa_current_subject.elem(previous_subjects));
      
      if (similarity_sum < epsilon){
        // if the similarities bewteen t-th subject and all previous subjects are 0
        arma::vec unique_previous_classes = unique(previous_classes);
        log_pmf += log(t-1) - log(m_0+t-1) - log(unique_previous_classes.size());
      }else{
        // if the similarities bewteen t-th subject and all previous subjects in the same class are 0
        // then add epsilon to avoid log_pmf_t being exactly -Inf
        if (similarity_sum_same_class < epsilon){ similarity_sum_same_class += epsilon; }
        // if both the similarity_sum and similarity_sum_same_class are non-zero
        log_pmf += log(t-1) - log(m_0+t-1) + log(similarity_sum_same_class) - log(similarity_sum);
      }
    }else{
      // if the class number of t-th subject is new 
      // i.e., it does not exist in previous t-1 subjects' class numbers
      log_pmf += log(m_0) - log(m_0+t-1);
    }  
  }
  
  // return log pmf of partition pi_n   
  return(log_pmf); 
}



// [[Rcpp::export]]
arma::vec update_beta_rcpp(const int k, const int q, arma::vec& index_class_k, 
                           arma::cube& gamma, arma::vec& e_kq, arma::mat& B_kq_inv, 
                           arma::cube& omega, const double sigma2, arma::cube& Y, 
                           arma::cube& X, arma::cube& H, arma::vec& J, 
                           const int Q, const int S, const int D){
  
  // update values of parameters beta_{kq}
  // 
  // args: k: the index of class number
  //       q: the index of depression item     
  //       index_class_k: index of subjects in class k
  //       gamma: kernel regression parameters 
  //       e_kq: mean of G_0 related to beta_{kq}
  //       B_kq_inv: inverse of covariance matrix of G_0 related to beta_{kq}
  //       omega: (n*J)*Q-dimensional normal correlation vector
  //       sigma2: variance of i.i.d normal error term
  //       data: Y, X, H, J, Q, S, D
  // returns:  beta_update_kq: new beta_{kq} after update
  
  arma::vec Xy_sum; Xy_sum.zeros(S); // S-dimension vector used for mean mu_n
  arma::mat X_sum; X_sum.zeros(S, S); // S*S matrix used for covariance matrix V_n
  arma::mat Y_q = Y.slice(q-1); 
  arma::mat omega_q = omega.slice(q-1);
  
  for (int iter=0; iter<index_class_k.size(); iter++){
    
    // index of i-th subject in class k 
    int i = index_class_k[iter]; 
    arma::vec Y_tilde = conv_to<vec>::from(Y_q.submat(i-1, 0, i-1, J[i-1]-1));
    
    arma::mat H_i = H.subcube(i-1, 0, 0, i-1, J[i-1]-1, D-1);
    arma::vec gamma_kq = gamma.subcube(k-1, q-1, 0, k-1, q-1, D-1);
    arma::vec H_gamma = conv_to<vec>::from(H_i * gamma_kq);
    
    arma::mat omega_iq = omega_q.submat(i-1, 0, i-1, J[i-1]-1);
    arma::vec omega_tilde = conv_to<vec>::from(omega_iq);
    
    // Y_tilde: J[i]-dimension vector used for calculate Xy_sum
    Y_tilde += - H_gamma - omega_tilde;
    
    if (J[i-1] == 1){
      arma::vec X_i = X.subcube(i-1, 0, 0, i-1, J[i-1]-1, S-1);
      Xy_sum += X_i * Y_tilde;
      X_sum += X_i * X_i.t();
    }else{
      arma::mat X_i = X.subcube(i-1, 0, 0, i-1, J[i-1]-1, S-1);
      Xy_sum += X_i.t() * Y_tilde;
      X_sum += X_i.t() * X_i;
    } 
  } 
  
  // Gaussian posterior distribution with mean mu_n and covariance V_n
  arma::mat V_n = inv(X_sum / sigma2 + B_kq_inv);
  arma::vec mu_n = V_n * (Xy_sum / sigma2 + B_kq_inv * e_kq);
  // Update beta_{kq}, k = 1,2,...,r_n, q = 1,2,...,Q
  arma::vec beta_update_kq = conv_to< vec >::from(rmvn_rcpp(1, mu_n, V_n));
  return(beta_update_kq);
}



// [[Rcpp::export]]
double logpost_sigma_omega_rcpp(arma::cube& omega, arma::mat& Sigma_omega, 
                                const int n, arma::vec& J, const int Q, const double sigma2){
  
  // log-posterior of correlation matrix Sigma_omega
  //
  // args: omega: normal correlation term
  //       Sigma_omega: correlation matrix of omega
  //       sigma2: variance of i.i.d normal error term
  //       data: n, J, Q
  // returns: logpr: log-posterior of Sigma_omega 
  
  double logll = 0; // log-likelihood of omega with correlation matrix Sigma_omega
  
  arma::rowvec omega_ij(Q);
  arma::vec z_ij(Q);
  arma::mat rooti = trans(inv(trimatu(chol(sigma2 * Sigma_omega))));
  double constants = sum(log(rooti.diag())) - (static_cast<double>(Q)/2.0)*log2pi;
  
  for (int i=0; i<n; i++){
    for (int j=0; j<J[i]; j++){
      omega_ij = omega.subcube(i, j, 0, i, j, Q-1);
      // calculate log dmvn of omega under Sigma_omega
      z_ij = rooti * trans(omega_ij);
      logll += constants - 0.5*sum(z_ij % z_ij);
    }
  }
  // LKJ prior of correlation matrix Sigma_omega
  double logpr = logll + log(det(Sigma_omega));
  return(logpr);
}



// [[Rcpp::export]]
arma::cube update_omega_rcpp(arma::cube& Y, arma::vec& pi_n, arma::cube& beta, 
                             arma::cube& gamma, arma::mat& Sigma_omega_inv, 
                             const double sigma2, const int n, arma::cube& X, 
                             arma::cube& H, arma::vec& J, const int Q, const int S, const int D){
  
  // update normal correlation term omega
  //
  // args: pi_n, beta, gamma: partition and parameters
  //       Sigma_omega_inv: inverse of Sigma_omega
  //       sigma2: variance of i.i.d normal error term
  //       data: n, Y, X, H, J, Q, S, D
  // returns: omega_update: new omega after update
  
  arma::cube omega_update(n, max(J), Q); omega_update.fill(0);
  // omega_update.fill(NA_REAL); 
  
  int class_index;
  arma::mat beta_i(Q, S);
  arma::mat gamma_i(Q, D);
  arma::colvec X_ij(S);
  arma::colvec Xbeta_ij(Q);
  arma::colvec H_ij(D);
  arma::colvec h_ij(Q);
  arma::colvec Y_ij(Q);
  arma::colvec Y_tilde(Q);
  arma::vec mu_n(Q);
  arma::mat V_n(Q, Q);
  arma::vec omega_update_ij(Q);
  arma::mat V_0 = 1/sigma2*eye<mat>(Q, Q);
  
  for (int i=0; i<n; i++){
    // index of class for i-th subject
    class_index = pi_n[i]-1; 
    beta_i = beta.subcube(class_index, 0, 0, class_index, Q-1, S-1);
    gamma_i = gamma.subcube(class_index, 0, 0, class_index, Q-1, D-1);
    for (int j=0; j<J[i]; j++){
      X_ij = X.subcube(i, j, 0, i, j, S-1);
      Xbeta_ij = conv_to< colvec >::from(beta_i * X_ij);
      H_ij = H.subcube(i, j, 0, i, j, D-1);
      h_ij = conv_to< colvec >::from(gamma_i * H_ij); 
      Y_ij = Y.subcube(i, j, 0, i, j, Q-1);
      Y_tilde = Y_ij - Xbeta_ij - h_ij;
      
      // Gaussian posterior distribution 
      V_n = inv(V_0 + Sigma_omega_inv);
      mu_n = conv_to< vec >::from(1/sigma2 * V_n * Y_tilde);
      omega_update_ij = conv_to< vec >::from(rmvn_rcpp(1, mu_n, V_n));
      omega_update.subcube(i, j, 0, i, j, Q-1) = omega_update_ij;
    }
  }
  return(omega_update);
}



// [[Rcpp::export]]
arma::mat update_sigma_omega_rcpp(arma::cube& omega, arma::mat& Sigma_omega, 
                                  const int n, arma::vec& J, const int Q, const double sigma2){
  
  // update covariance matrix of omega
  // 
  // args: omega: normal correlation term
  //       Sigma_omega: correlation matrix of omega
  //       sigma2: variance of i.i.d normal error term
  //       data: n, J, Q
  // returns: Sigma_omega_update: new Sigma_omega after update
  
  // double step = 0.025; // step size of proposal distribution
  double step = 0.01; // step size of proposal distribution
  double epsilon = 1e-300; // avoid numerical issue
  double rho, rho_new, lower, upper, ratio;
  arma::mat Sigma_omega_update = Sigma_omega;
  arma::mat Sigma_omega_new(Q, Q);
  
  for (int i=0; i<Q; i++){
    for (int j=0; j<Q; j++){
      if (i<j){
        // previous correlation rho, (i,j)-element of Sigma_omega  
        rho = Sigma_omega_update(i,j);
        lower = max(rho-step, -1+epsilon); 
        upper = min(rho+step, 1-epsilon);
        
        // propose new correlation rho from uniform distribution 
        rho_new = R::runif(lower, upper); // genereate new rho from uniform distribution
        Sigma_omega_new = Sigma_omega_update;
        Sigma_omega_new(i,j) = rho_new; Sigma_omega_new(j,i) = rho_new;  
        
        // check whether the new proposed matrix is positive definite
        if (!(Sigma_omega_new.is_sympd())) { continue; }
        ratio = logpost_sigma_omega_rcpp(omega, Sigma_omega_new, n, J, Q, sigma2) - 
          logpost_sigma_omega_rcpp(omega, Sigma_omega_update, n, J, Q, sigma2);
        
        // accept or reject new proposed rho
        if (log(R::runif(0,1)) < ratio){ Sigma_omega_update(i,j) = rho_new; Sigma_omega_update(j,i) = rho_new; }
      }
    }
  }
  return(Sigma_omega_update);
}



// [[Rcpp::export]]
double loglikelihood_rcpp(arma::uvec& Sk_minus,
                          arma::mat& e_tilde, arma::cube& B_tilde_inv_k, 
                          arma::mat& X_tilde_i, arma::mat& Y_tilde_i, 
                          arma::cube& X_tilde_sum, arma::cube& Xy_tilde_sum, 
                          const double sigma2, const int J_max, const int J_i, 
                          const int Q, const int S, const int D){
  
  // loglikelihood for existing clusters when update parition pi_n
  // 
  // args: Sk_minus: subjects in k-th cluster withour i-th subject
  //       e_tilde, B_tilde_inv_k: hyper-parameters for beta and gamma
  //       sigma2: variance of i.i.d normal error term
  //       data: X_tilde_i, Y_tilde_i, X_tilde_sum, Xy_tilde_sum, J_max, J_i, Q, S, D
  // returns: logll: loglikelihood
  
  double logll = 0;
  
  arma::mat B_tilde_inv_kq(S+D, S+D), V_kq(S+D, S+D), Sigma2_tilde(J_max, J_max);
  arma::vec Xy_tilde_sum_kq(S+D), e_tilde_q(S+D), mu_kq(S+D), mu_tilde(J_max);
  
  arma::mat X_tilde_sum_k = sum(X_tilde_sum.slices(Sk_minus), 2)/sigma2;
  arma::mat Xy_tilde_sum_k = sum(Xy_tilde_sum.slices(Sk_minus), 2)/sigma2; Xy_tilde_sum_k = Xy_tilde_sum_k.t();
  
  for (int q=0; q<Q; q++){
    
    // B_tilde_inv_kq = B_tilde_inv_k.row(q); 
    B_tilde_inv_kq = B_tilde_inv_k.slice(q);
    V_kq = inv_sympd(B_tilde_inv_kq + X_tilde_sum_k);
    
    e_tilde_q = e_tilde.col(q);
    Xy_tilde_sum_kq = Xy_tilde_sum_k.col(q);
    mu_kq = V_kq * (B_tilde_inv_kq * e_tilde_q + Xy_tilde_sum_kq);
    
    mu_tilde = X_tilde_i * mu_kq; 
    Sigma2_tilde = X_tilde_i * V_kq * X_tilde_i.t();
    
    for (int j=0; j<J_i; j++){
      logll += R::dnorm(Y_tilde_i(j, q), mu_tilde(j), sqrt(Sigma2_tilde(j, j) + sigma2), true);
    }
  }
  
  return(logll);
}



// [[Rcpp::export]]
double loglikelihood0_rcpp(arma::mat& e_tilde, arma::cube& B_tilde_inv_k, 
                           arma::mat& X_tilde_i, arma::mat& Y_tilde_i, 
                           const double sigma2, const int J_max, const int J_i, 
                           const int Q, const int S, const int D){
  
  // loglikelihood for the new cluster when update parition pi_n
  // 
  // args: e_tilde, B_tilde_inv_k: hyper-parameters for beta and gamma
  //       sigma2: variance of i.i.d normal error term
  //       data: X_tilde_i, Y_tilde_i, J_max, J_i, Q, S, D
  // returns: logll: loglikelihood
  
  double logll = 0;
  
  arma::mat B_tilde_inv_kq(S+D, S+D), B_tilde_kq(S+D, S+D), Sigma2_tilde(J_max, J_max);
  arma::vec e_tilde_q(S+D), mu_tilde(J_max);
  
  for (int q=0; q<Q; q++){
    
    B_tilde_inv_kq = B_tilde_inv_k.slice(q);
    B_tilde_kq = inv_sympd(B_tilde_inv_kq);
    
    e_tilde_q = e_tilde.col(q);
    mu_tilde = X_tilde_i * e_tilde_q;
    Sigma2_tilde = X_tilde_i * B_tilde_kq * X_tilde_i.t();
    
    for (int j=0; j<J_i; j++){
      logll += R::dnorm(Y_tilde_i(j, q), mu_tilde(j), sqrt(Sigma2_tilde(j, j) + sigma2), true);
    }
  }
  
  return(logll);
}



// [[Rcpp::export]]
arma::vec update_pi_n_rcpp(arma::vec& pi_n, const double m_0, arma::uvec& sigma, arma::mat& Kappa, 
                           arma::mat& e_tilde, arma::cube& B_tilde_inv, arma::mat& B_index,
                           arma::cube& X_tilde, arma::cube& Y_tilde, arma::cube& X_tilde_sum, arma::cube& Xy_tilde_sum, 
                           const int n, const int r_n_max, const double sigma2,  
                           arma::vec& J, const int Q, const int S, const int D){
  
  // update parition pi_n 
  // 
  // args: pi_n: current partition, m_0: mass parameter,  sigma: permutation,
  //       Kappa: ST kernel matrix, whose (i,j)-th element is the similarity between subject i and j 
  //       e_tilde, B_tilde_inv: hyper-parameters for beta and gamma
  //       sigma2: variance of i.i.d normal error term
  //       data: X_tilde, Y_tilde, X_tilde_sum, Xy_tilde_sum, B_index, n, r_n_max, J, Q, S, D
  // returns: pi_n_update: new partition after update
  
  arma::vec pi_n_update = pi_n;
  
  int r_n;
  double logll, log_pmf_pi_n;  
  arma::vec pi_n_new;
  arma::uvec Sk, Sk_minus;
  arma::cube B_tilde_inv_k(S+D, S+D, Q);
  arma::uvec B_index_k(Q);
  arma::mat X_tilde_i(max(J), S+D), Y_tilde_i(max(J), Q);
  
  // update the clustering membership of each subjects 
  for (int i=0; i<n; i++){
    
    pi_n_new = pi_n_update;
    r_n = max(pi_n_new);
    X_tilde_i = X_tilde.slice(i);
    Y_tilde_i = Y_tilde.slice(i);
    
    if (r_n >= r_n_max){
      Rcpp::NumericVector prob_classes(r_n);
      Rcpp::NumericVector log_prob_classes(r_n);
      
      // if moving i-th subject to exsiting class k = 1,2,...r_n
      for (int k=0; k<r_n; k++){
        
        pi_n_new[i] = k+1;
        Sk = find(pi_n_new == (k+1));
        Sk_minus = Sk.elem(find(Sk != i));
        
        if (Sk_minus.size() > 0){
          
          B_index_k = conv_to< arma::uvec >::from(B_index.col(k));
          B_tilde_inv_k = B_tilde_inv.slices(B_index_k);
          logll = loglikelihood_rcpp(Sk_minus, e_tilde, B_tilde_inv_k, X_tilde_i, Y_tilde_i, 
                                     X_tilde_sum, Xy_tilde_sum, sigma2, max(J), J[i], Q, S, D);
        }else{ 
          
          B_index_k = conv_to< arma::uvec >::from(B_index.col(r_n_max-1));
          B_tilde_inv_k = B_tilde_inv.slices(B_index_k);
          logll = loglikelihood0_rcpp(e_tilde, B_tilde_inv_k, X_tilde_i, Y_tilde_i, sigma2, max(J), J[i], Q, S, D);
        }
        
        log_pmf_pi_n = logpmf_pi_n_rcpp(pi_n_new, m_0, sigma, Kappa); 
        log_prob_classes[k] = logll + log_pmf_pi_n; 
      }
      
      // transform log probability to normalized probability and avoid overflow issue 
      log_prob_classes = log_prob_classes - max(log_prob_classes);
      prob_classes = exp(log_prob_classes) / sum(exp(log_prob_classes)); 
      
      // update the class number for i-th subject
      Rcpp::IntegerVector sample_set = seq(1, r_n);
      pi_n_update[i] = Rcpp::RcppArmadillo::sample(sample_set, 1, true, prob_classes)[0]; 
    }else{
      
      Rcpp::NumericVector prob_classes(r_n+1);
      Rcpp::NumericVector log_prob_classes(r_n+1);
      
      // if moving i-th subject to a new empty class, i.e., k = 0, index = r_n+1 
      pi_n_new[i] = r_n+1;
      B_index_k = conv_to< arma::uvec >::from(B_index.col(r_n_max-1));
      B_tilde_inv_k = B_tilde_inv.slices(B_index_k);
      logll = loglikelihood0_rcpp(e_tilde, B_tilde_inv_k, X_tilde_i, Y_tilde_i, sigma2, max(J), J[i], Q, S, D);
      log_pmf_pi_n = logpmf_pi_n_rcpp(pi_n_new, m_0, sigma, Kappa); 
      log_prob_classes[r_n] = logll + log_pmf_pi_n; 
      
      // if moving i-th subject to exsiting class k = 1,2,...r_n
      for (int k=0; k<r_n; k++){
        
        pi_n_new[i] = k+1;
        Sk = find(pi_n_new == (k+1));
        Sk_minus = Sk.elem(find(Sk != i));
        
        if (Sk_minus.size() > 0){
          
          B_index_k = conv_to< arma::uvec >::from(B_index.col(k));
          B_tilde_inv_k = B_tilde_inv.slices(B_index_k);
          logll = loglikelihood_rcpp(Sk_minus, e_tilde, B_tilde_inv_k, X_tilde_i, Y_tilde_i, 
                                     X_tilde_sum, Xy_tilde_sum, sigma2, max(J), J[i], Q, S, D);
        }else{ 
          
          B_index_k = conv_to< arma::uvec >::from(B_index.col(r_n_max-1));
          B_tilde_inv_k = B_tilde_inv.slices(B_index_k);
          logll = loglikelihood0_rcpp(e_tilde, B_tilde_inv_k, X_tilde_i, Y_tilde_i, sigma2, max(J), J[i], Q, S, D);
        }
        
        log_pmf_pi_n = logpmf_pi_n_rcpp(pi_n_new, m_0, sigma, Kappa); 
        log_prob_classes[k] = logll + log_pmf_pi_n; 
      }
      
      // transform log probability to normalized probability and avoid overflow issue 
      log_prob_classes = log_prob_classes - max(log_prob_classes);
      prob_classes = exp(log_prob_classes) / sum(exp(log_prob_classes)); 
      
      // update the class number for i-th subject
      Rcpp::IntegerVector sample_set = seq(1, r_n+1);
      pi_n_update[i] = Rcpp::RcppArmadillo::sample(sample_set, 1, true, prob_classes)[0]; 
    }
  }
  
  return(pi_n_update);
}