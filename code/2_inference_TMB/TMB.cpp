#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
    Type nll;           // Define variable that holds the return value (neg. log. lik)
    nll = Type(1.2);    // Assign value 1.2; a cast is needed.

    DATA_VECTOR(x);     // Vector x(0),x(1),...,x(n-1), where n is the length of x
    Type tmp = x(1);
    Type nll = tmp * tmp; 
    PARAMETER_VECTOR(u);          // Latent random variable 
    Type nll = Type(0);           // Return value
    nll -= dnorm(u(0),0,1,true)   // Distributional assignment: u(0) ~ N(0,1) 
}