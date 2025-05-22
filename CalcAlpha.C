template< class T > T square( const T& x ) { return x*x; }

void CalcAlpha()
{ 
  const double y_in = -0.0015;
  const double y_out = 0.0015;
  const double y_mean = (y_out+y_in)/2;
  const double y_mid = y_mean + 0.0005;
  const double y_cube = (square(y_in) + square( y_out ) + y_in*y_out)/3;
  
  const double numerator = y_out*(y_mean - y_mid) + (y_cube - y_mid*y_mean);
  const double denominator = y_in*(y_mean - y_mid) + (y_cube - y_mid*y_mean);
  
  const double alpha = 
    ( y_out*(y_mean - y_mid) - (y_mid*y_mean - y_cube) )/
    (  y_in*(y_mean - y_mid) - (y_mid*y_mean - y_cube) );
  
  std::cout << "CalcAlpha y_cube: " << y_cube << std::endl;
  std::cout << "CalcAlpha numerator: " << numerator << std::endl;
  std::cout << "CalcAlpha denominator: " << denominator << std::endl;
  std::cout << "CalcAlpha alpha: " << alpha << std::endl;
}
