#define SOR_HARD_STOP 100

#include <iostream>
#include "american_put.h"

double ** american_put(int N, int M, float r, double sigma, double K, double x_bar, double T, double Omega, double epsilon_stop)
{
    double ** result_table;
    result_table = new double *[M+1];
    for(int i =0; i<M+1; i++)
    {
        result_table[i] = new double[N];
        for(int j = 0; j<N; j++)
            result_table[i][j] = 0;
    }
    
    gsl_matrix * A = gsl_matrix_alloc(N,N);
    gsl_matrix_set_zero(A);
    for(int i = 0; i < N; i++)
        gsl_matrix_set(A, i, i, 1.0);
    
    gsl_vector * s_vec = gsl_vector_alloc(N);
    gsl_vector_set_all(s_vec, 1.0);
    
    projected_SOR(s_vec, A, Omega, epsilon_stop);
    
    gsl_matrix_free(A);
    gsl_vector_free(s_vec);
    return result_table;
};


void projected_SOR(gsl_vector * x, gsl_matrix * A, double Omega, double epsilon_stop)
{
    gsl_vector * x_half = gsl_vector_alloc(x->size);
    gsl_vector_set_all(x_half, 0.5);
    
    gsl_vector * x_epsilon = gsl_vector_alloc(x->size);
    
    bool converged = false;
    // Dummy Loop. modify at will!
    for(int k = 0; k<= SOR_HARD_STOP; k++)
    {
        gsl_vector_memcpy(x_epsilon, x);
        
        // doing some dummy algorithm
        gsl_vector_add(x, x_half);
        gsl_vector_scale (x, 0.5);
        
        // measuring how x has changed
        gsl_vector_sub(x_epsilon, x); // subtracts. x_epsilon - x stores in x_epsilon
        
        gsl_vector_mul(x_epsilon, x_epsilon); // make everything positive
        
        std::cout<<"max epsilon "<<gsl_vector_max(x_epsilon)<<std::endl;
        
        // Deciding to stop the interation or not.
        if(gsl_vector_max(x_epsilon) < epsilon_stop*epsilon_stop)
        {
            converged = true;
            break;
        }
    }
    
    if(!converged)
        std::cout<<"Warning: projected SOR DID NOT CONVERGE!!!"<<std::endl;
    
    gsl_vector_free(x_half);
    gsl_vector_free(x_epsilon);
    
};