template <class T>
double uni_slice(double x0, T &ll, double w, 
		 int m, double lower,  double upper, double gx0) 
{
  
  // Find the log density at the initial point, if not already known
  if(ISNAN(gx0)) {
    gx0 = ll.g(x0);
  }

  // Determine the slice level, in log terms
  double logy = gx0 - exp_rand();

  // Find the initial interval to sample from
  double u = R::runif(0.0, w);
  double L = x0 - u;
  double R = x0 + (w - u); // should guarantee that x0 is in [L,R], 
                           //even with roundoff
  
  // Expand the interval until its ends are outside the slice, 
  // or the limit on steps is reached.

  if(m == INT_MAX) {
    while(true) {
      if(L <= lower) break;
      if(ll.g(L) <= logy) break;
      L -= w;
    }
    while(true) {
      if(R >= upper) break;
      if(ll.g(R) <= logy) break;
      R += w;
    }
  }
  else if(m > 1) {
    int J = floor(R::runif(0.0, m));
    int K = (m-1) - J;

    while(J > 0) {
      if(L <= lower) break;
      if(ll.g(L) <= logy) break;
      L = L - w;
      J = J - 1;
    }
    
    while(K > 0) {
      if(R >= upper) break;
      if(ll.g(R) <= logy) break;
      R = R + w;
      K = K - 1;
    }
  }

  // Shrink interval to lower and upper bounds.

  if(L < lower) L = lower;
  if(R > upper) R = upper;

  // Sample from the interval, shrinking it on each rejection.

  double x1, gx1;

  while(true) {
    x1 = R::runif(L, R);
    gx1 = ll.g(x1);
    if(gx1 >= logy) break;
    if(x1 > x0) 
      R = x1;
    else
      L = x1;
  }

  return(x1);

}

template <class T>
double uni_slice1(double x0, T &ll, double w, 
		 int m, double lower,  double upper, double gx0) 
{
  
  // Find the log density at the initial point, if not already known
  if(ISNAN(gx0)) {
    gx0 = ll.g1(x0);
  }

  // Determine the slice level, in log terms
  double logy = gx0 - exp_rand();

  // Find the initial interval to sample from
  double u = R::runif(0.0, w);
  double L = x0 - u;
  double R = x0 + (w - u); // should guarantee that x0 is in [L,R], 
                           //even with roundoff
  
  // Expand the interval until its ends are outside the slice, 
  // or the limit on steps is reached.

  if(m == INT_MAX) {
    while(true) {
      if(L <= lower) break;
      if(ll.g1(L) <= logy) break;
      L -= w;
    }
    while(true) {
      if(R >= upper) break;
      if(ll.g1(R) <= logy) break;
      R += w;
    }
  }
  else if(m > 1) {
    int J = floor(R::runif(0.0, m));
    int K = (m-1) - J;

    while(J > 0) {
      if(L <= lower) break;
      if(ll.g1(L) <= logy) break;
      L = L - w;
      J = J - 1;
    }
    
    while(K > 0) {
      if(R >= upper) break;
      if(ll.g1(R) <= logy) break;
      R = R + w;
      K = K - 1;
    }
  }

  // Shrink interval to lower and upper bounds.

  if(L < lower) L = lower;
  if(R > upper) R = upper;

  // Sample from the interval, shrinking it on each rejection.

  double x1, gx1;

  while(true) {
    x1 = R::runif(L, R);
    gx1 = ll.g1(x1);
    if(gx1 >= logy) break;
    if(x1 > x0) 
      R = x1;
    else
      L = x1;
  }

  return(x1);

}

template <class T>
double uni_slice2(double x0, T &ll, double w, 
		 int m, double lower,  double upper, double gx0) 
{
  
  // Find the log density at the initial point, if not already known
  if(ISNAN(gx0)) {
    gx0 = ll.g2(x0);
  }

  // Determine the slice level, in log terms
  double logy = gx0 - exp_rand();

  // Find the initial interval to sample from
  double u = R::runif(0.0, w);
  double L = x0 - u;
  double R = x0 + (w - u); // should guarantee that x0 is in [L,R], 
                           //even with roundoff
  
  // Expand the interval until its ends are outside the slice, 
  // or the limit on steps is reached.

  if(m == INT_MAX) {
    while(true) {
      if(L <= lower) break;
      if(ll.g2(L) <= logy) break;
      L -= w;
    }
    while(true) {
      if(R >= upper) break;
      if(ll.g2(R) <= logy) break;
      R += w;
    }
  }
  else if(m > 1) {
    int J = floor(R::runif(0.0, m));
    int K = (m-1) - J;

    while(J > 0) {
      if(L <= lower) break;
      if(ll.g2(L) <= logy) break;
      L = L - w;
      J = J - 1;
    }
    
    while(K > 0) {
      if(R >= upper) break;
      if(ll.g2(R) <= logy) break;
      R = R + w;
      K = K - 1;
    }
  }

  // Shrink interval to lower and upper bounds.

  if(L < lower) L = lower;
  if(R > upper) R = upper;

  // Sample from the interval, shrinking it on each rejection.

  double x1, gx1;

  while(true) {
    x1 = R::runif(L, R);
    gx1 = ll.g2(x1);
    if(gx1 >= logy) break;
    if(x1 > x0) 
      R = x1;
    else
      L = x1;
  }

  return(x1);

}