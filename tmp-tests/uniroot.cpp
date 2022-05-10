// Code modified from https://github.com/wch/r-source/blob/79298c499218846d14500255efd622b5021c10ec/src/library/stats/src/zeroin.c

#include <Rcpp.h>
using namespace Rcpp;

double f(double x) {
  return(x - 1);
}


// [[Rcpp::export]]
double uniroot(double ax, double bx, double tol, double maxit) {

  double a = ax, b = bx, c = a;
  double fa = f(a), fb = f(b), fc = fa;

  if (fa == 0) return a;
  if (fb == 0) return b;

  while (maxit--) {  // Main iteration loop

    double prev_step = b - a;
    double tol_act, p, q, new_step;

    if (fabs(fc) < fabs(fb)) {
      // Swap data for b to be the best approximation
      a = b;    b = c;    c = a;
      fa = fb;  fb = fc;  fc= fa;
    }
    tol_act = 2 * DBL_EPSILON * fabs(b) + tol / 2;
    new_step = (c - b) / 2;

    if (fabs(new_step) <= tol_act || fb == 0) return b;

    /* Decide if the interpolation can be tried	*/
    if (fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb)) {
      double t1,cb,t2;
      cb = c-b;
      if( a==c ) {		/* If we have only two distinct	*/
    /* points linear interpolation	*/
    t1 = fb/fa;		/* can only be applied		*/
    p = cb*t1;
    q = 1 - t1;
      }
      else {			/* Quadric inverse interpolation*/

    q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
    p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1) );
    q = (q-1) * (t1-1) * (t2-1);
      }
      if( p>(double)0 )		/* p was calculated with the */
    q = -q;			/* opposite sign; make p positive */
    else			/* and assign possible minus to	*/
    p = -p;			/* q				*/

    if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
    && p < fabs(prev_step*q/2) )	/* and isn't too large	*/
    new_step = p/q;			/* it is accepted
     * If p/q is too large then the
     * bisection procedure can
     * reduce [b,c] range to more
     * extent */
    }

    if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
    if( new_step > (double)0 )	/* than tolerance		*/
    new_step = tol_act;
    else
      new_step = -tol_act;
    }
    a = b;	fa = fb;			/* Save the previous approx. */
    b += new_step;	fb = f(b);	/* Do step to a new approxim. */
    if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
      /* Adjust c for it to have a sign opposite to that of b */
      c = a;  fc = fa;
    }
  }

  return (fabs(fb) < fabs(fa)) ? b : a;  // failed, returns best solution for now
}
