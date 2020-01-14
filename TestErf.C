#include <TF1.h>
#include <TH1.h>

#include <iostream>
#include <cmath>

Double_t FitFunction( Double_t* x, Double_t* par )
{
  const Double_t t0 = (x[0]-par[1]+par[3]*std::sqrt(12)/2)/par[2];
  const Double_t t1 = (x[0]-par[1]-par[3]*std::sqrt(12)/2)/par[2];
  return par[0]*( std::erf(t0) + std::erf(-t1) );
}

void TestErf( void )
{

  TF1* f = new TF1( "fitFunction", FitFunction, -10, 10, 4 );
  f->SetParameters( 1, 0, 0, 5 );
  // f->Draw();
  auto h = f->GetHistogram();
  h->Draw();
  std::cout << "RMS: " << h->GetRMS() << std::endl;

}
