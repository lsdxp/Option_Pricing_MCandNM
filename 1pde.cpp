#include <iostream>
#include <fstream>
#include <cmath>

#include "do_nm.h"

int main(void)
{
    std::ofstream output_file;
    output_file.open("/Users/shengdongliu/Desktop/out.csv");
    int N = 20;
    int M = 10;
    double r = 0.02;
    double sigma = .2;
    double K = 20;
    double S_max = 100;
    double T = 20;
    
//    int N = 50;
//    int M = 20;
//    double r = 0.1;
//    double sigma = .4;
//    double K = 10;
//    double S_max = 40;
//    double T =0.25 ;
    
    
    //double ** do_euler_explicit(int N, int M, float r, double sigma, double K, double S_max, double T)
    double ** result_table = NULL;
    //result_table = do_euler_explicit(N, M, r, sigma, K, S_max, T);
    
    //result_table = do_euler_implicitet(N, M, r, sigma, K, S_max, T);
    
    result_table = do_crank_nicolson(N, M, r, sigma, K, S_max, T);
    
    
    std::cout<<"The Result!\n";
    for(int i =0; i < M; i++)
    {
        for(int j = 0; j<N-1; j++)
        {
            std::cout<<result_table[i][j]<<" ";
            output_file<<result_table[i][j]<<',';
        }
        std::cout<<result_table[i][N-1]<<std::endl;
        output_file<<result_table[i][N-1]<<"\n";
    }
    
}