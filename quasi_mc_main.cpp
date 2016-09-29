#include <cstdlib> // *** Thanks to Leonhard Gruenschloss and Mike Giles   ***
#include <cmath>   // *** for pointing out the change in new g++ compilers ***
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#define pi 3.1415926
clock_t t;

using namespace std;

// ----- SOBOL POINTS GENERATOR BASED ON GRAYCODE ORDER -----------------
// INPUT:
//   N         number of points  (cannot be greater than 2^32)
//   D         dimension  (make sure that the data file contains enough data!!)
//   dir_file  the input file containing direction numbers
//
// OUTPUT:
//   A 2-dimensional array POINTS, where
//
//     POINTS[i][j] = the jth component of the ith point,
//
//   with i indexed from 0 to N-1 and j indexed from 0 to D-1
//
// ----------------------------------------------------------------------

double get_uniform()
{
    return (((double) random())/(pow(2.0, 31.0)-1.0));
}

double **sobol_points(unsigned N, unsigned D, char *dir_file)
{
    ifstream infile(dir_file,ios::in);
    if (!infile) {
        cout << "Input file containing direction numbers cannot be found!\n";
        exit(1);
    }
    char buffer[1000];
    infile.getline(buffer,1000,'\n');
    
    // L = max number of bits needed
    unsigned L = (unsigned)ceil(log((double)N)/log(2.0));
    
    // C[i] = index from the right of the first zero bit of i
    unsigned *C = new unsigned [N];
    C[0] = 1;
    for (unsigned i=1;i<=N-1;i++) {
        C[i] = 1;
        unsigned value = i;
        while (value & 1) {
            value >>= 1;
            C[i]++;
        }
    }
    
    // POINTS[i][j] = the jth component of the ith point
    //                with i indexed from 0 to N-1 and j indexed from 0 to D-1
    double **POINTS = new double * [N];
    for (unsigned i=0;i<=N-1;i++) POINTS[i] = new double [D];
    for (unsigned j=0;j<=D-1;j++) POINTS[0][j] = 0;
    
    // ----- Compute the first dimension -----
    
    // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
    unsigned *V = new unsigned [L+1];
    for (unsigned i=1;i<=L;i++) V[i] = 1 << (32-i); // all m's = 1
    
    // Evalulate X[0] to X[N-1], scaled by pow(2,32)
    unsigned *X = new unsigned [N];
    X[0] = 0;
    for (unsigned i=1;i<=N-1;i++) {
        X[i] = X[i-1] ^ V[C[i-1]];
        POINTS[i][0] = (double)X[i]/pow(2.0,32); // *** the actual points
        //        ^ 0 for first dimension
    }
    
    // Clean up
    delete [] V;
    delete [] X;
    
    
    // ----- Compute the remaining dimensions -----
    for (unsigned j=1;j<=D-1;j++) {
        
        // Read in parameters from file
        unsigned d, s;
        unsigned a;
        infile >> d >> s >> a;
        unsigned *m = new unsigned [s+1];
        for (unsigned i=1;i<=s;i++) infile >> m[i];
        
        // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
        unsigned *V = new unsigned [L+1];
        if (L <= s) {
            for (unsigned i=1;i<=L;i++) V[i] = m[i] << (32-i);
        }
        else {
            for (unsigned i=1;i<=s;i++) V[i] = m[i] << (32-i);
            for (unsigned i=s+1;i<=L;i++) {
                V[i] = V[i-s] ^ (V[i-s] >> s);
                for (unsigned k=1;k<=s-1;k++)
                    V[i] ^= (((a >> (s-1-k)) & 1) * V[i-k]);
            }
        }
        
        // Evalulate X[0] to X[N-1], scaled by pow(2,32)
        unsigned *X = new unsigned [N];
        X[0] = 0;
        for (unsigned i=1;i<=N-1;i++) {
            X[i] = X[i-1] ^ V[C[i-1]];
            POINTS[i][j] = (double)X[i]/pow(2.0,32); // *** the actual points
            //        ^ j for dimension (j+1)
        }
        
        // Clean up
        delete [] m;
        delete [] V;
        delete [] X;
    }
    delete [] C;
    
    return POINTS;
}




class twovalue{
    double z1;
    double z2;
public:
    void set_values (double,double);
    double get1() { return z1;}
    double get2() { return z2; }
};
void twovalue::set_values (double x, double y) {
    z1= x;
    
    z2 = y;
}

twovalue get_gaussian(double sobol1,double sobol2)
{
    double R,V,Z1,Z2;
    twovalue G;
    R=-2.0*log(sobol1);
    V=2*pi*sobol2;
    Z1=(sqrt(R)*cos(V));
    Z2=(sqrt(R)*sin(V));
    G.set_values(Z1, Z2);
    return G ;
}




int main(int argc, char **argv)
{
    if (argc != 4) {
        cout << endl << "input format: sobol N D FILENAME" << endl << endl;
        cout << "The program prints the first N sobol points in D dimensions." << endl;
        cout << "The points are generated in graycode order." << endl;
        cout << "The primitive polynomials and initial direction numbers are" << endl
        << "given by the input file FILENAME." << endl << endl;
        return 0;
    }
    
//    int N = atoi(argv[1]);
//    int D = atoi(argv[2]);
    
//T = 2, r = 0.05, σ = 0.5, S = K = 2, q = 0.
    
    double Si,r,q,sigma,T,K,St,temp,xbar,ybar,m,delta;
//    r=10.0/100;
//    sigma=0.2;
//    Si=100;
//    T=1.0;
//    K=100;
//    m=50;
//    q=0.0;
    
    r=5.0/100;
    sigma=0.5;
    Si=2;
    T=2.0;
    K=2;
    m=70;
    q=0.0;

    delta=T/m;
    int l=20;

    
    double callprice,sum;
    temp=pow(sigma,2.0);
    
    double interval_upperbound, interval_lowerbound,quantile,estimated_standard_error,average_S,price_sum,width_of_confidence_interval,S0,Splus1,Splus2;
    quantile=1.96;
    double n=5000.0;
    int N=int(n);
    int D=int(m);
    
    double **P = sobol_points(N,D,argv[3]);
    
    // display points
    //cout << setprecision(20);
    //cout << setiosflags(ios::scientific) << setprecision(10);
//    for (unsigned i=0;i<=N-1;i++) {
//        for (unsigned j=0;j<=D-1;j++) cout << P[i][j] << " " ;
//        cout << endl;
//    }
//    cout << endl;
    

    double average[l];
    
    double shift[D];
    
    
    
    
    double g1,g2;
    twovalue G;
    

    
    
    
    
//
    for (int ll=0;ll<20;ll++)
    {
     sum=0.0;
        for (int i=0;i<D;i++)
            shift[i]=get_uniform();
        
        for (unsigned i=0;i<=N-1;i++)
        {
            for (unsigned j=0;j<=D-1;j++)
                P[i][j]=(P[i][j]+shift[j])-floor(P[i][j]+shift[j]);
        }
    
    
    for(int i=0;i<n;i++)
    {
        S0=Si;
        price_sum=0.0;
        for(int j=0;j<D;)
        {
            G=get_gaussian(P[i][j],P[i][j+1]);
            g1=G.get1();
            g2=G.get2();
            
            Splus1=S0*exp( (r-0.5*sigma*sigma) *delta +sigma*sqrt(delta)*g1);
            Splus2=Splus1*exp((r-0.5*sigma*sigma) *delta +sigma*sqrt(delta)*g2);
            price_sum=price_sum+Splus1+Splus2;
            S0=Splus2;
            j=j+2;
        }
        average_S=price_sum/m;
        callprice=exp(-r*T)*max((average_S-K),0.0);
        sum+=callprice;
 
    }
    
    average[ll]=sum/(n);
    }
    int batch=20;
    double total=0.0;
    double k=2.0;
    double avr=0.0;
    width_of_confidence_interval=1.0;
    
    
    
    for(int i=0;i<batch;i++)
    {
        total+=average[i];
        if(i==0)
        {
            xbar=average[i];
            ybar=pow(average[i], 2.0);
            continue;
        }
        else
        {
            xbar=(1-(1/k))*xbar+(1/k)*average[i];
            ybar=(1-(1/k))*ybar+(1/k)*pow(average[i],2.0);
            k=k+1.0;
        }
    }
//    double var=0.0;
//    for (int i=0;i<20;i++)
//        
//        var+=(average[i]-7.16176)*(average[i]-7.16176);
//    var=var/l;
//    
//    cout<<"Asian option call price : "<<average[0]<<endl;
//    cout<<sqrt(var)<<endl;
    
    estimated_standard_error=sqrt((1.0/(batch-1))*(ybar-pow(xbar,2.0)));
    
    avr=total/double(batch);
    //estimated_standard_error=sigma/sqrt(n);
    interval_lowerbound=avr-quantile*estimated_standard_error;
    interval_upperbound=avr+quantile*estimated_standard_error;
    width_of_confidence_interval=interval_upperbound-interval_lowerbound;
    
    
    
    cout<<"number of trials :"<<n*batch<<std::setprecision(6)<<endl;
    cout<<"dimension m :"<<m<<endl;
    cout<<"Asian option call price : "<<avr<<endl;
    cout<<"estimated standard error :"<<estimated_standard_error<<endl;
    cout<<"95% confidence interval : ["<<interval_lowerbound<<" , "<<interval_upperbound<<"]"<<endl;
    cout<<"width of confidence interval :"<<width_of_confidence_interval<<endl;
    t = clock() - t;
    cout<<"time consumed :"<<(float)(t)/float(CLOCKS_PER_SEC)<<" seconds"<<endl;
    cout<<"efficiency :"<<pow(estimated_standard_error,2.0)*(float)(t)/float(CLOCKS_PER_SEC)<<endl;
    
  
}