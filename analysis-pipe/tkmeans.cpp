// constructing priority queues
#include <iostream>       
#include <queue>         
#include <vector>         
#include <functional>   
#include <tuple>
// maths
#include <math.h> 
#include <RcppArmadillo.h>


using namespace Rcpp;
const double log2pi = std::log(2.0 * M_PI);

typedef std::tuple<double, int, int>  di_pair;

class CompareDist
{
public:
  bool operator()(di_pair n1, di_pair n2) {
    return std::get<0>(n1)<std::get<0>(n2);
  }
};


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov) {
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);    
}

// [[Rcpp::export]]
arma::vec dmvnorm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool log = false) { 
  arma::vec distval = Mahalanobis(x,  mean, sigma);
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
  
  if (log) { 
    return(logretval);
  } else { 
    return(exp(logretval));
  }
}


// [[Rcpp::export]]
arma::vec distCentre(int n_cluster, arma::mat ui, arma::mat centres, double lambda, int d){
  //set up sigma
  arma::vec sigmaDiag(d);  sigmaDiag.fill(lambda);
  arma::mat sigma = diagmat(sigmaDiag);
  
   //store distances
  arma::vec dist =  arma::zeros(n_cluster);
  
  for (int j=0; j<n_cluster; j++)
  {
    dist.at(j) = 1.0/n_cluster * dmvnorm_arma(ui, centres.row(j), sigma, false)[0];
  }
  //Rcpp::Rcout << "size: " << size(dist) << std::endl;
  return(dist);
}

// [[Rcpp::export]]
arma::vec distCentre2(int n_cluster, arma::mat ui, arma::mat centres, double lambda, int d){
  
  //store distances
  arma::vec dist =  arma::zeros(n_cluster);
  
  for (int j=0; j<n_cluster; j++)
  {
    dist.at(j) = sqrt(arma::sum(arma::pow(centres.row(j) - ui, 2.0)));
  }
  
  
  return(dist);
}

// [[Rcpp::export]]
arma::vec distCentre3(int n_cluster, arma::mat ui, arma::mat centres, double lambda, int d){
  arma::vec dist = arma::sqrt(arma::sum(arma::pow(centres - repmat(ui,n_cluster,1), 2.0),1));
  return(dist);
}


// [[Rcpp::export]]
int max_index(arma::vec x){
  
  //returns index of biggest value in x.
  Rcpp::NumericVector y = NumericVector(x.begin(),x.end());
  Rcpp::NumericVector::iterator it =  std::max_element(y.begin(), y.end());
  
  //Rcpp::Rcout << "it  " << it <<  " y " << y.begin() <<  " diff " << std::distance(y.begin(), it) << std::endl;
  return (std::distance(y.begin(), it));
}


// [[Rcpp::export]]
int min_index(arma::vec x){
  
  //returns index of biggest value in x.
  Rcpp::NumericVector y = NumericVector(x.begin(),x.end());
  Rcpp::NumericVector::iterator it =  std::min_element(y.begin(), y.end());
  
  //Rcpp::Rcout << "it  " << it <<  " y " << y.begin() <<  " diff " << std::distance(y.begin(), it) << std::endl;
  return (std::distance(y.begin(), it));
}


// [[Rcpp::export]]
arma::mat init_centres(arma::mat M, int n_cluster){
  int n = size(M)[0];
  int d = size(M)[1];
  
  //randomly pick centres
  arma::mat B = arma::randi<arma::mat>(n_cluster,1,arma::distr_param(0,n));
  
  arma::mat cent =  arma::zeros(n_cluster,d);
  
  for(int i=0; i<n_cluster; i++){
    //Rcpp::Rcout << "i = " << i << " n = " << B.at(i,0) <<  std::endl;
    
    cent.row(i) = M.row(B.at(i,0));
  }
  
  return cent;
}


// [[Rcpp::export]]
arma::mat tkmeans_lowmem(arma::mat& M, int n_cluster , double alpha, int nstart, int iter){
  int n = size(M)[0];
  int d = size(M)[1];
  
  Rcpp::Rcout << "n = " << n << " d = " << d <<  " k = " << n_cluster << std::endl;
  Rcpp::Rcout << "n starts = " << nstart << " max iterations " << iter <<  " alpha = " << alpha << std::endl;
  
  arma::mat centres = init_centres(M, n_cluster); 
  arma::mat means = centres;
  
  double lambda = 1;
  double tol = 0.001;
  double diff = 1;
  int m=0;
  
  // cluster membership record
  //std::vector<int> cluster_membership;
  //cluster_membership.reserve(n); 
  
  //keep list of smallest using priority queue
  std::priority_queue< di_pair, std::vector<di_pair>, CompareDist  >  smallest;
  
  int queue_len;
  if(alpha != 0.0){
    queue_len = floor((n-1)*(1-alpha))+1;
  }else{
    queue_len = 0;
  }
  
  
  //Rcpp::Rcout << "number of outliers = " <<  queue_len   << std::endl;
  //Rcpp::Rcout << "Begin main loop..."  << std::endl;
  
  
  while(m<iter & diff > tol){
    
    //Rcpp::Rcout << m << " of " << iter  << " iterations" << std::endl;
    
    arma::mat new_centres = arma::zeros(n_cluster, d);
    arma::mat centre_members = arma::zeros(n_cluster, 1);
    
    for(int i=0; i<n; i++)
    {
        arma::mat Dc = distCentre2(n_cluster, M.row(i), centres, lambda, d);
        //distCentre(n_cluster, M.row(i), centres, lambda, d).print();
        
        int cluster = min_index(Dc);
        
        di_pair temp = di_pair(Dc.at(cluster), cluster, i); 
        
        //new_centres.row(cluster) =  new_centres.row(cluster)+ M.row(i);
        //centre_members(cluster,0) =  centre_members(cluster,0) + 1;
        
        if(queue_len > 0){
          if(smallest.size() < queue_len){
             smallest.push(temp);
          }else{
            if(std::get<0>(smallest.top()) > std::get<0>(temp)){
              //Rcpp::Rcout << "out: "   <<  std::get<0>(smallest.top()) <<  " in: " << std::get<0>(temp) << std::endl;
              smallest.pop();
              smallest.push(temp);
            }
          }
        }
        //Rcpp::Rcout << "size of heap:"  <<  smallest.size() << std::endl;
    }
    
   
  //remove all the outliers
    if(queue_len > 0){
        while (!smallest.empty())
        {
          new_centres.row(std::get<1>(smallest.top())) = new_centres.row(std::get<1>(smallest.top())) + M.row(std::get<2>(smallest.top()));
          
          centre_members(std::get<1>(smallest.top()),0) = centre_members(std::get<1>(smallest.top()),0) + 1;
          
          //Rcpp::Rcout << " " << std::get<0>(smallest.top()) <<  " " << std::get<1>(smallest.top()) << " " << std::get<2>(smallest.top()) << std::endl;
          
          smallest.pop();
        }
    }
    
    if(centre_members.min() > 0) {
        means = new_centres.each_col() / centre_members;
    }else{
       stop("Empty cluster");
    }
    
    //centre_members.print();
    //Rcpp::Rcout << diff << std::endl;
    
    diff = arma::accu(arma::abs(centres-means));
    centres = means;
    m++;
  }
  
  return means;
}


// [[Rcpp::export]]
arma::mat scale_lowmem(arma::mat& M){
  arma::rowvec means = mean(M);
  arma::rowvec sds = stddev(M);
  
  for(int i=0; i< M.n_cols; i++){
    M.col(i) -= means.at(i);
    M.col(i) /= sds.at(i);
  }
  
  arma::mat means_sds(2,means.n_elem);
  means_sds.row(0) = means;
  means_sds.row(1) = sds;
  
  return means_sds;
  }


// [[Rcpp::export]]
double  BIC_lowmem(arma::mat& centres, arma::mat& data){
          int x = data.n_rows;
          int m = centres.n_cols;
          int k = centres.n_rows;
          double PI_val = log(1./k);
          
          //Rcpp::Rcout <<  x<< ' ' << m <<' ' << k << std::endl;
          
          
          arma::mat temp_row(1, k);
          
          
          double log_like_accum = 0.;
          
          for(int i=0;i<x;i++){
            
            for(int j=0;j<k;j++){

              temp_row.at(0,j) = 
                
                -0.5*(arma::accu(arma::pow(data.row(i),2.)) 
                + arma::accu(arma::pow(centres.row(j),2.))  
                - 2.*dot(data.row(i),centres.row(j)))
                - (m/2.)*log(2.*M_PI) - log(m)/2. + PI_val;
  
            }
            
            double temp_max = temp_row.max();
            
            
            temp_row = temp_row - temp_max;
            
            double log_like = temp_max + log(arma::accu(arma::exp(temp_row)));
            //Rcpp::Rcout <<  log_like << std::endl;
            log_like_accum += log_like;
          }
          
          //Rcpp::Rcout <<  log_like_accum << std::endl;
          
          
  return -2*log_like_accum  + log(x)*(m*k + k - 1);
}

// [[Rcpp::export]]
arma::uvec tmeansClust_lowmem(arma::mat& data, arma::mat& centres){
  int n_cluster =  centres.n_rows;
  double lambda = 1.0;
  int d = size(data)[0];
  
  arma::uvec clusters = arma::zeros<arma::uvec>(d);
  
  for(int i=0;i<d;i++){
    arma::mat Dc = distCentre2(n_cluster, data.row(i), centres, lambda, d);
    int cluster = min_index(Dc);
    clusters.at(i) = cluster+1;//going back to 1 indexed R
  }
  
  return clusters;
}

