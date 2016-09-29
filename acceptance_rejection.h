//
//  f.h
//  hw3.2
//
//  Created by Shengdong  Liu on 2/27/16.
//  Copyright © 2016 Shengdong  Liu. All rights reserved.
//

#ifndef f_h
#define f_h
#define c 0.918938533204672
using namespace std;

double get_vvalue(int j)
{
    double v[15];
    v[0]=1.253314137315500;
    v[1]=0.6556795424187985;
    v[2]=0.4213692292880545;
    v[3]=0.3045902987101033;
    v[4]=0.2366523829135607;
    v[5]=0.1928081047153158;
    v[6]=0.1623776608968675;
    v[7]=0.1401041834530502;
    v[8]=0.1231319632579329;
    v[9]=0.1097872825783083;
    v[10]=0.09902859647173193;
    v[11]=0.09017567550106468;
    v[12]=0.08276628650136917;
    v[13]=0.0764757610162485;
    v[14]=0.07106958053885211;
    return (v[j]);
}


double get_uniform()
{
    return (((double) random())/(pow(2.0, 31.0)-1.0));
}

double min(double a, int b) {
    return (b > a )? a:b;
}



double max(double a, double b) {
    return (b < a )? a:b;
}




#endif /* f_h */
