#include <iostream>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>


double max(double a, double b) {
    return (b < a )? a:b;
}

double ** do_euler_explicit(int N, int M, float r, double sigma, double K, double S_max, double T)
{
//--------
//    /**************************************************************************************************************/
//    /***** Below is an example, replace with whatever you want to use to build this method ************************/
//    
//    double delta_t=T/M;
//    double delta_S=S_max/N;
//    double middle1=-(1+delta_t*r)-delta_t*sigma*sigma/(delta_S*delta_S);
////    double middle2=;
//    double plus1=-((sigma*sigma*0.5*(1/(2*delta_S)-(1/(delta_S*delta_S)))-(r/(2*delta_S)))*delta_t);
////    double plus2=0.5*sigma*sigma*delta_t;
//    double minus1=-((sigma*sigma*0.5*(-1/(2*delta_S)-(1/(delta_S*delta_S)))+(r/(2*delta_S)))*delta_t);
////    double minus2=0.5*delta_t*sigma*sigma;
//    
//    double Si[N];
//    for(int i =0;i<N;i++)
//        Si[i]=(i+1)*delta_S;
//    
//    
//    double ** result_table;
//    result_table = new double *[M];
//    for(int i =0; i<M; i++)
//    {
//        result_table[i] = new double[N];
//        for(int j = 0; j<N; j++)
//        result_table[i][j] = 0;
//    }
//    
//    // Build Example Matrix
//    gsl_matrix * t_m = gsl_matrix_alloc(N,N);
//    gsl_matrix_set_zero(t_m);
//    
//    
//    for(int i = 0; i < N; i++)
//    gsl_matrix_set(t_m, i, i,middle1);
//    
//    for(int i = 0; i < N-1; i++)
//    gsl_matrix_set(t_m, i+1, i, minus1);
//    
//    for(int i = 0; i <N-1; i++)
//    gsl_matrix_set(t_m, i, i+1, plus1);
//    
//    for(int i = 0; i < N; i++)
//    {
//        for(int j = 0; j < N; j++)
//        std::cout<<gsl_matrix_get(t_m, i, j)<<" ";
//        std::cout<<std::endl;
//    }
//    
//    // Solution Vector
//    gsl_vector * s_v_old = gsl_vector_alloc(N);
//    gsl_vector * s_v_new = gsl_vector_alloc(N);
//    
//    
//    gsl_vector_set_zero(s_v_old);
//    gsl_vector_set_zero(s_v_new);
//   
//    
//    for(int i = 0; i < N; i++)
//        gsl_vector_set(s_v_old,i,max((Si[i]-K),0.0));
//    
//    
//    //example initial value
//    
//   //    gsl_vector_set(s_v_old, N/2, middle1+middle2*gsl_vector_get(call, N/2));
////    
////    
////    gsl_vector_set(s_v_old, N/2 - 1, minus1*gsl_vector_get(call, -1+N/2)+minus2*gsl_vector_get(call, -1+N/2)*gsl_vector_get(call,-1+N/2));
////    
////    gsl_vector_set(s_v_old, N/2 + 1, plus1*gsl_vector_get(call, 1+N/2)+plus2*gsl_vector_get(call, 1+N/2)*gsl_vector_get(call, 1+N/2));
////
//    
//    std::cout<<"New\n";
//    for(int i = 0; i < N; i++)
//    std::cout<<gsl_vector_get(s_v_new, i)<<" ";
//    std::cout<<std::endl;
//    
//    std::cout<<"Old\n";
//    for(int i = 0; i < N; i++)
//    std::cout<<gsl_vector_get(s_v_old, i)<<" ";
//    std::cout<<std::endl;
//    
//    // printing out the initial condition
//    std::cout<<"Starting\n";
//    for(int i = 0; i < N; i++)
//    std::cout<<gsl_vector_get(s_v_old, i)<<" ";
//    std::cout<<std::endl;
//    
//    // this is the loop that steps through time
//    for(int t_j = 0; t_j < M; t_j++)
//    {
//        // This is what actually does the matrix vector multiplication.
//        gsl_blas_dgemv(CblasNoTrans, 1.0, t_m, s_v_old, 0.0, s_v_new);
//        
//        for(int i = 0; i < N; i++)
//        { if(gsl_vector_get(s_v_new, i)>(S_max-K))
//                gsl_vector_set(s_v_new, i,(S_max-K));
//            if(gsl_vector_get(s_v_new, i)<0.0)
//                gsl_vector_set(s_v_new, i,0.0);
//        }
////            if(gsl_vector_get(s_v_new, N-1)>(S_max-K))
////                gsl_vector_set(s_v_new, N-1,(S_max-K));
////            if(gsl_vector_get(s_v_new, 0)<0.0)
////                gsl_vector_set(s_v_new, 0,0.0);
//    
//        
//        for(int i = 0; i < N; i++)
//        std::cout<<gsl_vector_get(s_v_new, i)<<" ";
//        std::cout<<std::endl;
//        
//        
//        memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
//        
//        gsl_vector_swap(s_v_new, s_v_old);
//    }
//    
//    gsl_matrix_free(t_m);
//    gsl_vector_free(s_v_old);
//    gsl_vector_free(s_v_new);
//    
//    /***** Above is an example, replace with whatever you want to use to build this method ************************/
//    /**************************************************************************************************************/
//    
//    return result_table;
//------------
    
    
    double ** result_table;
    result_table = new double *[M];
    for(int i =0; i<M; i++)
    {
        result_table[i] = new double[N];
        for(int j = 0; j<N; j++)
            result_table[i][j] = 0;
    }
    
    double delta_t=T/M;
    double delta_S=2*log(S_max)/N;
    double middle1=1-delta_t*r-delta_t*sigma*sigma/(delta_S*delta_S);
    double a1,a2,a3;
    a1=sigma*sigma*0.5;
    a2=1/(2*delta_S);
    a3=1/(delta_S*delta_S);
    
   
    
    double plus1=(a1*(a3-a2)+r*a2)*delta_t;
    
    double minus1=(a1*(a3+a2)-r*a2)*delta_t;

    
    double Si[N];
        for(int i =0;i<N;i++)
            Si[i]=-log(S_max)+(i+1)*delta_S;
    // Build Example Matrix
    gsl_matrix * t_m = gsl_matrix_alloc(N,N);
    gsl_matrix_set_zero(t_m);
    for(int i = 0; i < N; i++)
        gsl_matrix_set(t_m, i, i, middle1);

    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m, i+1, i, minus1);
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m, i, i+1, plus1);
    
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            std::cout<<gsl_matrix_get(t_m, i, j)<<" ";
        std::cout<<std::endl;
    }
    
    
    
    // Solution Vector
    gsl_vector * s_v_old = gsl_vector_alloc(N);
    gsl_vector * s_v_new = gsl_vector_alloc(N);
    
    gsl_vector_set_zero(s_v_old);
    gsl_vector_set_zero(s_v_new);
    
    //example initial value
    
//    gsl_vector_set(s_v_old, N/2, 2.5);
//    gsl_vector_set(s_v_old, N/2 - 1, 1.5);
//    gsl_vector_set(s_v_old, N/2 + 1, 1.5);
    
    for(int i = 0; i < N; i++)
        gsl_vector_set(s_v_old,i,max((exp(Si[i])-K*exp(-r*T)),0.0));
    
    std::cout<<"New\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_new, i)<<" ";
    std::cout<<std::endl;
    
    std::cout<<"Old\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    // printing out the initial condition
    std::cout<<"Starting\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    // this is the loop that steps through time
    for(int t_j = 0; t_j < M; t_j++)
    {
        // This is what actually does the matrix vector multiplication.
        gsl_blas_dgemv(CblasNoTrans, 1.0, t_m, s_v_old, 0.0, s_v_new);
        gsl_vector_set(s_v_new,N-1,(S_max-K*exp(-r*(T-(t_j+1)*delta_t))));
        gsl_vector_set(s_v_new,0,0.0);
//        for(int i=0;i<N;i++)
//        {  if((gsl_vector_get(s_v_new, i)!=0.0)&(gsl_vector_get(s_v_new, i)<0.0))
//        {
//            gsl_vector_set(s_v_new, i,0.0);
//            break;
//        }
//            if((gsl_vector_get(s_v_new, i)!=0.0)&(gsl_vector_get(s_v_new, i)>(S_max-K)))
//            { gsl_vector_set(s_v_new, i,S_max-K);
//                break;
//            }
//            
//            
//        }
//        for(int i=N-1;i>0;i--)
//              {
//                if((gsl_vector_get(s_v_new, i)!=0.0)&(gsl_vector_get(s_v_new, i)>(S_max-K*exp(-r*(i+1)*delta_t))))
//                { gsl_vector_set(s_v_new, i,S_max-K*exp(-r*(i+1)*delta_t));
//                    break;
//                }
//                  if((gsl_vector_get(s_v_new, i)!=0.0)&(gsl_vector_get(s_v_new, i)<0.0))
//                  {
//                      gsl_vector_set(s_v_new, i,0.0);
//                      break;
//                  }
//              }
    
        for(int i = 0; i < N; i++)
            std::cout<<gsl_vector_get(s_v_new, i)<<" ";
        std::cout<<std::endl;
        
        memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
        
        gsl_vector_swap(s_v_new, s_v_old);
    }
    
    gsl_matrix_free(t_m);
    gsl_vector_free(s_v_old);
    gsl_vector_free(s_v_new);
    
    /***** Above is an example, replace with whatever you want to use to build this method ************************/
    /**************************************************************************************************************/
    
    return result_table;
    
    
};

double ** do_euler_implicitet(int N, int M, float r, double sigma, double K, double S_max, double T)
{
    /**************************************************************************************************************/
    /***** Below is an example, replace with whatever you want to use to build this method ************************/
    

    double delta_t=T/M;
    double delta_S=2*log(S_max)/N;
    double middle1=1+delta_t*r+delta_t*sigma*sigma/(delta_S*delta_S);
    double a1,a2,a3;
    a1=sigma*sigma*0.5;
    a2=1/(2*delta_S);
    a3=1/(delta_S*delta_S);
    
    
    
    double plus1=(a1*(-a3+a2)-r*a2)*delta_t;
    
    double minus1=(a1*(-a3-a2)+r*a2)*delta_t;
    
    
    double Si[N];
    for(int i =0;i<N;i++)
        Si[i]=-log(S_max)+(i+1)*delta_S;

    
    double ** result_table;
    result_table = new double *[M];
    for(int i =0; i<M; i++)
    {
        result_table[i] = new double[N];
        for(int j = 0; j<N; j++)
        result_table[i][j] = 0;
    }
    
    // Build Example Matrix
    gsl_matrix * t_m = gsl_matrix_alloc(N,N);
    gsl_matrix_set_zero(t_m);
    gsl_matrix_set_zero(t_m);
    for(int i = 0; i < N; i++)
        gsl_matrix_set(t_m, i, i, middle1);
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m, i+1, i, minus1);
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m, i, i+1, plus1);
    
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            std::cout<<gsl_matrix_get(t_m, i, j)<<" ";
        std::cout<<std::endl;
    }
    
    // Solution Vector
    gsl_vector * s_v_old = gsl_vector_alloc(N);
    gsl_vector * s_v_new = gsl_vector_alloc(N);
    
    gsl_vector_set_zero(s_v_old);
    gsl_vector_set_zero(s_v_new);
    
    for(int i = 0; i < N; i++)
        gsl_vector_set(s_v_old,i,max((exp(Si[i])-K*exp(-r*T)),0.0));

    
    //example initial value
//    gsl_vector_set(s_v_old, N/2, -(middle1+middle2*2*delta_S/N));
//    gsl_vector_set(s_v_old, N/2 - 1, -(minus1*2*delta_S/N+minus2*2*delta_S/N*2*delta_S/N));
//    gsl_vector_set(s_v_old, N/2 + 1, -(plus1*2*delta_S/N+plus2*2*delta_S/N*2*delta_S/N));
//

    std::cout<<"New\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_new, i)<<" ";
    std::cout<<std::endl;
    
    std::cout<<"Old\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    // printing out the initial condition
    std::cout<<"Starting\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    // this is the loop that steps through time
    for(int t_j = 0; t_j < M; t_j++)
    {
        // This is what actually does the matrix vector multiplication.
//        gsl_blas_dgemv(CblasNoTrans, 1.0, t_m, s_v_old, 0.0, s_v_new);
        gsl_permutation * p = gsl_permutation_alloc (N);
        int s;
        gsl_linalg_LU_decomp (t_m, p, &s);
        gsl_linalg_LU_solve (t_m, p, s_v_old,s_v_new);
        
        for(int i = 0; i < N; i++)
            std::cout<<gsl_vector_get(s_v_new, i)<<" ";
        std::cout<<std::endl;
        
        gsl_vector_set(s_v_new,N-1,(S_max-K*exp(-r*(T-(t_j+1)*delta_t))));
        gsl_vector_set(s_v_new,0,0.0);
        
        
        
        memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
        
        gsl_vector_swap(s_v_new, s_v_old);
    }
    
    gsl_matrix_free(t_m);
    gsl_vector_free(s_v_old);
    gsl_vector_free(s_v_new);
    
    /***** Above is an example, replace with whatever you want to use to build this method ************************/
    /**************************************************************************************************************/
    
    return result_table;
};


double ** do_crank_nicolson(int N, int M, float r, double sigma, double K, double S_max, double T)
{
    /**************************************************************************************************************/
    /***** Below is an example, replace with whatever you want to use to build this method ************************/
    
    double delta_t=T/M;
    double delta_S=2*log(S_max)/N;
    double middle1=1+0.5*delta_t*r+0.5*delta_t*sigma*sigma/(delta_S*delta_S);
    double a1,a2,a3;
    a1=sigma*sigma*0.5;
    a2=1/(2*delta_S);
    a3=1/(delta_S*delta_S);
    
    double plus1=0.5*((a1*(-a3+a2)-r*a2)*delta_t);
    
    double minus1=0.5*((a1*(-a3-a2)+r*a2)*delta_t);
    
    
    double middle1_right=1-0.5*delta_t*r-0.5*delta_t*sigma*sigma/(delta_S*delta_S);
    
    double plus1_right=0.5*((a1*(a3-a2)+r*a2)*delta_t);
    
    double minus1_right=0.5*((a1*(a3+a2)-r*a2)*delta_t);
    
    
    
    double Si[N];
    for(int i =0;i<N;i++)
    
        Si[i]=-log(S_max)+(i+1)*delta_S;
    
    
    
    double ** result_table;
    result_table = new double *[M];
    for(int i =0; i<M; i++)
    {
        result_table[i] = new double[N];
        for(int j = 0; j<N; j++)
            result_table[i][j] = 0;
    }
    
    
    
    // Build Example Matrix
    gsl_matrix * t_m = gsl_matrix_alloc(N,N);
    gsl_matrix_set_zero(t_m);
    
    for(int i = 0; i < N; i++)
        gsl_matrix_set(t_m, i, i, middle1);
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m, i+1, i, minus1);
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m, i, i+1, plus1);

    
   
    gsl_matrix * c_m = gsl_matrix_alloc(N,N);
    gsl_matrix_set_zero(c_m);
    
    gsl_matrix * t_m_right = gsl_matrix_alloc(N,N);
    gsl_matrix_set_zero(t_m_right);
    
    for(int i = 0; i < N; i++)
        gsl_matrix_set(t_m_right, i, i, middle1_right);
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m_right, i+1, i, minus1_right);
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m_right, i, i+1, plus1_right);
    
    
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            std::cout<<gsl_matrix_get(t_m, i, j)<<" ";
        std::cout<<std::endl;
    }
    
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            std::cout<<gsl_matrix_get(t_m_right, i, j)<<" ";
        std::cout<<std::endl;
    }
    
    // Solution Vector
    gsl_vector * s_v_old = gsl_vector_alloc(N);
    gsl_vector * s_v_new = gsl_vector_alloc(N);
    gsl_vector * s_v_temp = gsl_vector_alloc(N);

    
    gsl_vector_set_zero(s_v_old);
    gsl_vector_set_zero(s_v_new);
    gsl_vector_set_zero(s_v_temp);

    
   
    //example initial value
//    gsl_vector_set(s_v_old, N/2, -0.5*(middle1_left+middle2_left*2*delta_S/N));
//    gsl_vector_set(s_v_old, N/2 - 1, -0.5*(minus1_left*2*delta_S/N+minus2_left*2*delta_S/N*2*delta_S/N));
//    gsl_vector_set(s_v_old, N/2 + 1, -0.5*(plus1_left*2*delta_S/N+plus2_left*2*delta_S/N*2*delta_S/N));
//
    
    for(int i = 0; i < N; i++)
        
    {
        gsl_vector_set(s_v_old,i,max((exp(Si[i])-K*exp(-r*T)),0.0));
    }
    
    std::cout<<"New\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_new, i)<<" ";
    std::cout<<std::endl;
    
    std::cout<<"Old\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    // printing out the initial condition
    std::cout<<"Starting\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    
    
    gsl_permutation * p = gsl_permutation_alloc (N);
    int s;
    gsl_linalg_LU_decomp (t_m, p, &s);
    
    // this is the loop that steps through time
    for(int t_j = 0; t_j < M; t_j++)
    {
        // This is what actually does the matrix vector multiplication.
        //        gsl_blas_dgemv(CblasNoTrans, 1.0, t_m, s_v_old, 0.0, s_v_new);
        //gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,t_m_right, s_v_old,0.0,c_m);
        
        
        //gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, t_m_right, s_v_old,0.0, c_m);
        //  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0, t_m_right,s_v_old,0.0, s_v_old);
        
        gsl_blas_dgemv (CblasNoTrans, 1.0, t_m_right, s_v_old, 0.0, s_v_temp);
        
        //s_v_old=gsl_matrix_mul_matrix_by_matrix ( t_m_right ,s_v_old);
        
        gsl_linalg_LU_solve (t_m, p, s_v_temp,s_v_new);
        
        gsl_vector_set(s_v_new,N-1,(S_max-K*exp(-r*(T-(t_j+1)*delta_t))));
        gsl_vector_set(s_v_new,0,0.0);
        
        for(int i = 0; i < N; i++)
            std::cout<<gsl_vector_get(s_v_new, i)<<" ";
        std::cout<<std::endl;
        
        memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
        
        gsl_vector_swap(s_v_new, s_v_old);
        
    }
    
    gsl_matrix_free(t_m);
    gsl_vector_free(s_v_old);
    gsl_vector_free(s_v_new);

    
    
    /***** Above is an example, replace with whatever you want to use to build this method ************************/
    /**************************************************************************************************************/
    return result_table;
};