/*
see http://www.ibiblio.org/e-notes/Splines/Basis.htm
to understand the following.

for n parameters to be fitted 
for instance for phase, n_degree_phase,
what is the bspline step (s) for a given range of validity [t_min,t_max] ?
for order k, you need k parameters (points d'ancrage) to model the value at t.
 n*s = t_max-t_min+(k-1)*s
 s=(t_max-t_min)/(n-k+1)
*/


static double sqr(double x) {return x*x;}


static double Bspline3(double t, int i) {
  if( (t<i) || (t>(i+3)) ) return 0; 
  if(t<(i+1)) return 0.5*sqr(t-i);
  if(t<(i+2)) return 0.5*( (i+2-t)*(t-i)+(t-i-1)*(i+3-t) );
  return 0.5*sqr(i+3-t);
}

// recursive evaluation of the spline
// see http://www.ibiblio.org/e-notes/Splines/Basis.htm
static double Bspline(int k, double t, int i) {
  switch(k){
  case 3:
    {
      if( (t<i) || (t>(i+3)) )
	return 0;
      if(t<(i+1))
	return 0.5*sqr(t-i);
      if(t<(i+2))
	return 0.5*( (i+2-t)*(t-i)+(t-i-1)*(i+3-t) );
      return 0.5*sqr(i+3-t);
    }
    break;
  case 1:
    {
      if( (t>=i) && (t<(i+1)) )
	return 1.0;
      else
	return 0.0;
    }
    break;
  case 2:
    {
      if( (t<i) || (t>(i+2)) )
	return 0;
      if(t<(i+1))
	return t-i;
      return i+2-t;
    }
    break;
    
  default:
   return (Bspline(k-1,t,i)*(t-i)+Bspline(k-1,t,i+1)*(i+k-t))/(k-1.0);
   break;
  }
}
