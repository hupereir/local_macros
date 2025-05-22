#include <TCanvas.h>
#include <TF1.h>

#include <cmath>

//_____________________________________________________
double square( double x ) { return x*x; }

//_____________________________________________________
double Gaus( double* x, double* par )
{                         
  const double xx = x[0];
  const double amplitude = par[0];
  const double mean = par[1];
  const double sigma = par[2];
  return std::exp( -square(xx/sigma)/2 )/(sigma*std::sqrt(M_PI*2));
}

//_____________________________________________________
double integral( double xmin, double xmax, double sigma )
{ return ( std::erf( xmax/(std::sqrt(2)*sigma ) ) - std::erf( xmin/(std::sqrt(2)*sigma ) ) )/2; }

//_____________________________________________________
void TestGaus( void )
{
  
  double mean = 0;
  double sigma = 1;
  
  auto f = new TF1( "Gaus", Gaus, -5, 5, 3 );
  f->SetParameter(0,1);
  f->SetParameter(1,mean);
  f->SetParameter(2,sigma);
  
  auto cv = new TCanvas( "cv", "cv", 900, 900 );
  f->Draw();
  
  // full instegral
  std::cout << "Integral: " << f->Integral( -5, 5 ) << std::endl;
  
  // partial integral
  const double a = 0;
  const double b = 2;
  std::cout << "Integral: " << f->Integral( a, b ) << " calc: " << integral( a, b, sigma ) << std::endl;
  
}
