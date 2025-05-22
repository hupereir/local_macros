template< class T > T square( const T& x ) { return x*x; }

void CalcAlpha_new()
{
  const double y_in = -0.0015;
  const double y_out = 0.0015;
  
  const double r_in = 8e-4;
  const double r_out = 25e-4;
  
  const int nseg = 100;
  
  double sumyr2 = 0;
  double sumr2 = 0;
  
  for( int i = 0; i < nseg; ++i )
  {
    
    // calculate center
    const double y = y_in + (0.5 + i)*(y_out - y_in)/nseg; 
    const double r = r_in + (y-y_in)*(r_out - r_in)/(y_out-y_in);
    sumyr2 += y*square(r);
    sumr2 += square(r);
  }

  const double y_mean = sumyr2/sumr2;
  std::cout << "CalcAlpha_new - y_mean: " << y_mean << std::endl;
  
}
