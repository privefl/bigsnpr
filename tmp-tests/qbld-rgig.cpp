#include <Rcpp.h>

inline double square(double x) { return x * x; }

inline double mode(double lambda, double omega) {
  double one_minus_lambda = 1 - lambda;
  return omega / (one_minus_lambda + ::sqrt(square(one_minus_lambda) + square(omega)));
}

double rgig_noshift(double lambda, int check, double omega, double alpha)
{
  double xm,nc,ym,um,s,t,U,V,X;

  t = 0.5*(lambda-1);
  s = 0.25*omega;

  xm = mode(lambda,omega);
  nc = t*log(xm) - s*(xm + 1/xm);
  ym = ((lambda+1) + sqrt((lambda+1)*(lambda+1) + omega*omega))/omega;
  um = exp(0.5*(lambda+1)*log(ym) - s*(ym + 1/ym) - nc);

  do{
    U = um*::unif_rand();
    V = ::unif_rand();
    X = U/V;
  } while ((log(V)) > (t*log(X) - s*(X+1/X)- nc));

  return (check==1) ? (alpha/X) : (alpha*X);
}

double rgig_shift(double lambda, int check, double omega, double alpha)
{
  double xm,nc,s,t,U,V,X,a,b,c,p,q,fi,fak,y1,y2,uplus,uminus;

  t = 0.5*(lambda-1);
  s = 0.25*omega;

  xm = mode(lambda,omega);
  nc = t*log(xm) - s*(xm + 1/xm);

  a = -(2*(lambda+1)/omega +xm);
  b = (2*(lambda-1)*xm/omega -1);
  c = xm;

  p = b - a*a/3;
  q = (2*a*a*a)/27 - (a*b)/3 + c;

  fi = acos(-q/(2*sqrt(-p*p*p/27)));
  fak = 2*sqrt(-p/3);
  y1 = fak*cos(fi/3) - a/3;
  y2 = fak*cos(fi/3 + (4./3.)*M_PI) - a/3;

  uplus = (y1-xm)*exp(t*log(y1) - s*(y1 + 1/y1) -nc);
  uminus =  (y2-xm)*exp(t*log(y2) - s*(y2 + 1/y2) -nc);

  do{
    U = uminus + ::unif_rand() * (uplus-uminus);
    V = ::unif_rand();
    X = U/V + xm;
  } while ((X<=0) || (log(V) > (t*log(X) - s*(X+1/X)- nc)));

  return (check==1) ? (alpha/X) : (alpha*X);
}


double rgig_conc(double lambda, int check, double omega, double alpha)
{
  double A1, A2, A3, Atot,k0,k1,k2,xm,x0,a,U,V,X,hx;

  if(lambda >=1 || omega > 1)
    Rcpp::stop("Invalid parameters: lambda or omega");

  xm = mode(lambda,omega);
  x0 = omega/(1-lambda);

  k0 = exp((lambda-1)*log(xm) - 0.5*omega*(xm + 1/xm));
  A1 = k0*x0;

  if (x0 >= 2/omega) {
    k1 = 0;
    A2 = 0;
    k2 = pow(x0,lambda-1);
    A3 = k2*2*exp(-omega*x0/2)/omega;
  } else {
    k1 = exp(-omega);
    A2 = (lambda==0) ? (k1*log(2/(omega*omega))) :
      ((k1/lambda)*(pow(2/omega,lambda) - pow(x0,lambda)));
    k2 = pow(2/omega,lambda-1);
    A3 = k2*2*exp(-1.0)/omega;
  }

  Atot = A1 + A2 + A3;

  do {

    V = Atot * ::unif_rand();

    if (V <= A1) {
      X = x0 * V / A1;
      hx = k0;
    } else {

      V -= A1;

      if (V <= A2) {
        if (lambda == 0) {
          X = omega * exp(exp(omega)*V);
          hx = k1 / X;
        } else {
          X = pow(pow(x0, lambda) + (lambda / k1 * V), 1/lambda);
          hx = k1 * pow(X, lambda-1);
        }
      } else {
        V -= A2;
        a = (x0 > 2/omega) ? x0 : 2/omega;
        X = -2/omega * log(exp(-omega/2 * a) - omega/(2*k2) * V);
        hx = k2 * exp(-omega/2 * X);
      }
    }

    U = hx * ::unif_rand();

  } while(log(U) > (lambda-1)*log(X) - omega/2 * (X+1/X));

  return (check==1) ? (alpha/X) : (alpha*X);
}


// [[Rcpp::export]]
double rgig_one(double lambda, double chi, double psi) {

  // expecting chi > 0 and psi > 0
  // in my application, lambda is always the same

  int check = 0;

  if (lambda < 0) {
    lambda = -lambda;
    check = 1;
  }

  double alpha = ::sqrt(psi / chi);
  double omega = ::sqrt(psi * chi);

  if (lambda > 2 || omega > 3) {
    return rgig_shift(lambda, check, omega, alpha);          // RoU shift
  } else if (lambda >= (1 - 2.25 * chi * psi) || omega > 0.2) {
    return rgig_noshift(lambda, check, omega, alpha);        // RoU no shift
  } else {
    return rgig_conc(lambda, check, omega, alpha);           // log-concave
  }
}

/*** R
hist(replicate(1e5, rgig_one(lambda = 0.5, chi = 0.1, psi = 30)))
hist(qbld::rgig(1e5, lambda = 0.5, a = 0.1, b = 30))
*/
