
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>

//_____________________________________________________
double square( double x ) { return x*x; }

//_____________________________________________________
double Gaus( double* x, double* par )
{
  const double amplitude = par[0];
  const double mean = par[1];
  const double sigma = par[2];

  const double xx = (x[0]-mean);
  return amplitude*std::exp( -square(xx/sigma)/2 );
}

//_________________________________________________________
double SquareStrip( double* x, double* par )
{
  const double amplitude = par[0];
  const double mean = par[1];
  const double sigma = par[2];

  const double xx = (x[0]-mean);
  if( std::abs( xx ) < sigma/2 ) return amplitude;
  else return 0;
}

//_________________________________________________________
double TriangularStrip( double* x, double* par )
{
  const double amplitude = par[0];
  const double mean = par[1];
  const double sigma = par[2];

  const double xx = (x[0]-mean);
  if( std::abs( xx ) < sigma ) return amplitude*(sigma-std::abs(xx))/sigma;
  else return 0;
}

//_________________________________________________________
void DrawSharing()
{
 
  set_style( false );
  
  // dummy histogram
  auto h = new TH1F( "h", "",  100, -3, 3 );
  h->SetMinimum(0);
  h->SetMaximum(1.5);
  
  auto cv = new TCanvas( "cv", "cv", 900, 700 );
  h->Draw();
  
  // strips
  const double pitch = 1;
  const double width = 0.98;
  for( int i = 0; i < 3; ++i )
  {
    double mean = -1 + i*pitch;
    const TString fname = Form( "strip_%i", i );
    auto f = new TF1( fname, SquareStrip, mean-2*pitch, mean+2*pitch, 3 );
    // auto f = new TF1( fname, TriangularStrip, mean-2*pitch, mean+2*pitch, 3 );
    f->SetNpx( 1000 );
    f->SetParameter(0,1);
    f->SetParameter(1,mean);
    f->SetParameter(2,width);
    f->SetLineColor( i == 1 ? 1:17);
    f->Draw("same");
  }
 
  // draw gaussian
  auto f = new TF1( "gaus", Gaus, -5, 5, 3 );
  f->SetNpx( 1000 );
  f->SetParameter(0,1.3);
  f->SetParameter(1,1.1);
  f->SetParameter(2,0.4);
  f->SetLineColor(2);
  
  f->Draw("same");
  
  cv->SaveAs( "Figures/ChargeSharing_square.pdf" );
  // cv->SaveAs( "Figures/ChargeSharing_triangular.pdf" );
    
  
}
