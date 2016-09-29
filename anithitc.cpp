#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;
#define pi 3.141592626535

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


clock_t t;

double get_uniform()
{
    return (((double) random())/(pow(2.0, 31.0)-1.0));
}



twovalue get_gaussian()

{
    double R,V,Z1,Z2;
    twovalue G;
    R=-2.0*log(get_uniform());
    V=2*pi*get_uniform();
    Z1=(sqrt(R)*cos(V));
    Z2=(sqrt(R)*sin(V));
    G.set_values(Z1, Z2);
    
    
    return G ;
}

double max(double a, double b) {
    return (b < a )? a:b;
}

int main (int argc, char* argv[])
{
    double Si,r,q,sigma,T,K,St,Stminus,temp,average,xbar,ybar;
    r=0.3866/100;
    sigma=29.79/100;
    Si=1868.99;
    T=1.0/52.0;
    K=1870;
    q=2.32/100;
    
    double callprice,sum,callprice1,callprice2;
    // sum=0.0;
    
    temp=pow(sigma,2.0);
    
    double interval_upperbound, interval_lowerbound,quantile,estimated_standard_error,width_of_confidence_interval;
    quantile=1.96;
    double n=200000000.0;
    //width_of_confidence_interval=1.0;
    
    sum=0.0;
    double k=2.0;
    double g1,g2;
    twovalue G;
    for(int i=1;i<=n;i++)
    {
        G=get_gaussian();
        g1=G.get1();
        g2=G.get2();
        
        St=Si*exp ( (r-q-(temp/2.0) )*T+sigma*sqrt(T)*g1);
        Stminus=Si*exp((r-q-(temp/2.0) )*T-sigma*sqrt(T)*g2);
        callprice1=exp(-r*T)*(max((St-K),0.0));
        callprice2=exp(-r*T)*(max((Stminus-K),0.0));
        callprice=0.5*(callprice1+callprice2);

        
        sum+=callprice;
        if(i==1)
        {
            xbar=callprice;
            ybar=pow(callprice, 2.0);
            continue;
        }
        
        xbar=(1-(1/k))*xbar+(1/k)*callprice;
        ybar=(1-(1/k))*ybar+(1/k)*pow(callprice,2.0);
        k=k+1.0;
        
    }
    estimated_standard_error=sqrt((1/(n-1))*(ybar-pow(xbar,2.0)));
    average=sum/(n);
    
    // estimated_standard_error=sigma/sqrt(n);
    interval_lowerbound=average-quantile*estimated_standard_error;
    interval_upperbound=average+quantile*estimated_standard_error;
    width_of_confidence_interval=interval_upperbound-interval_lowerbound;
    
    
    
    cout<<"number of trials :"<<n<<endl;
    cout<<"European option call price : "<<average<<endl;
    cout<<"estimated standard error :"<<estimated_standard_error<<endl;
    cout<<"95% confidence interval : ["<<interval_lowerbound<<" , "<<interval_upperbound<<"]"<<endl;
    cout<<"width of confidence interval :"<<width_of_confidence_interval<<endl;
    t = clock() - t;
    cout<<"time consumed :"<<(float)(t)/float(CLOCKS_PER_SEC)<<" seconds"<<endl;
    cout<<"efficiency :"<<pow(estimated_standard_error,2.0)*(float)(t)/float(CLOCKS_PER_SEC)<<endl;
    return 0;
    
}
