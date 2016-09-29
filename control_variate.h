//
//  h1.h
//  MC_hw3
//
//  Created by Shengdong  Liu on 2/26/16.
//  Copyright © 2016 Shengdong  Liu. All rights reserved.
//

#ifndef h1_h
#define h1_h

#define m1 2147483647
#define m2 2145483479
#define a12 63308
#define a13 -183326
#define a21 86098
#define a23 -539608
#define q12 33921
#define q13 11714
#define q21 24919
#define q23 3976
#define r12 12979
#define r13 2883
#define r21 7417
#define r23 2071
#define Invmp1 4.656612873077393e-10

#include <iostream>
#include <cmath>

//int x10, x11, x12, x20, x21, x22;
int x10=283646;
int x11=214748;
int x12=2146;
int x20=214548;
int x21=483478;
int x22=478;
int Random()
{
    int h, p12, p13, p21, p23;
    // std::cout<<Z<<std::endl;
    //    x10=((double) random())/(pow(2.0, 31.0)-1.0)*2147483645+1;
    //    x11=((double) random())/(pow(2.0, 31.0)-1.0)*2147483645+1;
    //    x12=((double) random())/(pow(2.0, 31.0)-1.0)*2147483645+1;
    //    x20=((double) random())/(pow(2.0, 31.0)-1.0)*2145483477+1;
    //    x21=((double) random())/(pow(2.0, 31.0)-1.0)*2145483477+1;
    //    x22=((double) random())/(pow(2.0, 31.0)-1.0)*2145483477+1;
    
    h = x10/q13;
    // std::cout<<random()<<std::endl;
    
    
    p13 = -a13*(x10-h*q13)-h*r13;
    h = x11/q12;
    p12 = a12*(x11-h*q12)-h*r12;
    if(p13<0)
        p13 = p13+m1;
    if(p12<0)
        p12 = p12+m1;
    x10 = x11;
    x11 = x12;
    x12 = p12-p13;
    if(x12<0)
        x12 = x12+m1;
    h = x20/q23; p23 = -a23*(x20-h*q23)-h*r23;
    h = x22/q21; p21 = a21*(x22-h*q21)-h*r21;
    if(p23<0)
        p23 = p23+m2;
    if(p21<0)
        p21 = p21+m2;
    x20 =x21;
    x21 =x22;
    x22 =p21 - p23;
    if(x22 <0)
        x22 =x22 +m2;
    if (x12<x22)
        return (x12-x22+m1);
    else
        return (x12-x22);
}
double Uniform01()
{
    int Z;
    Z=Random();
    //std::cout<<Z<<std::endl;
    if(Z==0)
        Z=m1;
    return(Z*Invmp1);
}

double get_uniform()
{
    return (((double) random())/(pow(2.0, 31.0)-1.0));
}


#endif /* h1_h */
