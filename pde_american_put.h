#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

double ** american_put(int N, int M, float r, double sigma, double K, double x_bar, double T, double Omega, double epsilon_stop);
void projected_SOR(gsl_vector * x, gsl_matrix * A, double Omega, double epsilon_stop);