//
//  do_nm.h
//  pde0404
//
//  Created by Shengdong  Liu on 4/4/16.
//  Copyright © 2016 Shengdong  Liu. All rights reserved.
//

#ifndef do_nm_h
#define do_nm_h


double ** do_euler_explicit(int N, int M, float r, double sigma, double K, double S_max, double T);
double ** do_euler_implicitet(int N, int M, float r, double sigma, double K, double S_max, double T);
double ** do_crank_nicolson(int N, int M, float r, double sigma, double K, double S_max, double T);


#endif /* do_nm_h */
