#ifndef UTIL_H
#define UTIL_H

#include <complex>

extern const double PI;

std::complex<double> ** new_cmat (int m, int n); 
std::complex<double> ** czeros (int m, int n); 
std::complex<double> ** eyec (int m, int n); 
std::complex<double> * copyv (const std::complex<double> *orig, int n); 
void copyv(const std::complex<double> *orig, std::complex<double> *vnew, int n);
std::complex<double> ** copym (std::complex<double> **orig, int m, int n); 
std::complex<double> ** copym (std::complex<double> **orig, int n); 
void copym (std::complex<double> **orig, std::complex<double> **mnew, int m, 
            int n); 
void copym (std::complex<double> **orig, std::complex<double> **mnew, int n); 
std::complex<double> ** vec_to_mat(const std::complex<double> *x, int m, int n);
void deletem (std::complex<double> **A, int m); 

std::complex<double> ** conj_transpose (std::complex<double> **A, int m, int n);
std::complex<double> ** conj_transpose (std::complex<double> **A, int n);

std::complex<double> ** sum (std::complex<double> **A, std::complex<double> **B,
                             int n); 
void sum_ip (std::complex<double> **A, std::complex<double> **B, int n); 
std::complex<double> ** diff(std::complex<double> **A, std::complex<double> **B,
                             int n); 
void diff_ip (std::complex<double> **A, std::complex<double> **B, int n); 

std::complex<double> inner_prod (const std::complex<double> *x, 
                                 const std::complex<double> *y, int n); 
std::complex<double> ** outer_prod (const std::complex<double> *x, 
                                    const std::complex<double> *y, int n); 
std::complex<double> * prod (std::complex<double> **A, 
                             const std::complex<double> *x, int m, int n); 
std::complex<double> * prod (std::complex<double> **A, 
                             const std::complex<double> *x, int n); 
std::complex<double> ** prod (std::complex<double> **A, 
                              std::complex<double> **B, int m, int n, int p); 
std::complex<double> ** prod (std::complex<double> **A, 
                              std::complex<double> **B, int n); 
void prod_ip (std::complex<double> **A, std::complex<double> *x, int n); 
void prod_ip (std::complex<double> **A, std::complex<double> **B, int m, int n);
void prod_ip (std::complex<double> **A, std::complex<double> **B, int n);
std::complex<double> ** kronecker(std::complex<double> **A,
                                  std::complex<double> **B,
                                  int m, int n, int p, int q);
std::complex<double> ** mpow (std::complex<double> **A, int e, int n); 
void mpow_ip (std::complex<double> **A, int e, int n); 

bool is_unitary (std::complex<double> **U, int n);

int gcd (int a, int b); 
bool is_prime (int n); 
int int_root (int n); 

double runif ();
int rdunif (int a, int b);

std::string int_to_bin (int n, int len); 
int get_bit (int n, int index); 

void normalize (std::complex<double> *coeffs, int n); 


#endif 







