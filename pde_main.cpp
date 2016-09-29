#include <iostream>
#include <fstream>

#include "american_put.h"

int main(void)
{
    std::ofstream output_file;
    output_file.open("/Users/shengdongliu/Desktop/IE525MC/pde_americanput/out.csv");
    int N = 4;
    int M = 10;
    double r = 0.02;
    double sigma = .2;
    double K = 20;
    double S_max = 50;
    double T = 20;
    double Omega = 1.0;
    double epsilon_stop = 0.001;
    double ** result_table = NULL;
    // double ** american_put(int N, int M, float r, double sigma, double K, double x_bar, double T, double Omega, double epsilon_stop);
    result_table = american_put(N, M, r, sigma, K, S_max, T, Omega, epsilon_stop);
    std::cout<<"The Result!\n";
    for(int j = 0; j<N-1; j++)
    {
        std::cout<<j<<" ";
        output_file<<j<<',';
    }
    std::cout<<N-1<<std::endl;
    output_file<<N-1<<"\n";
    
    for(int i =0; i < M + 1; i++)
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