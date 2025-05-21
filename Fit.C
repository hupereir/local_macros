#ifndef _Fit_C_
#define _Fit_C_

#include <TF1.h>
#include <TH1.h>
#include <TMath.h>

#include <cmath>

//_______________________________________
template<class T> T square( T x ) { return x*x; }

//_______________________________________
class Result
{

  public:

  Result()
  {}

  Result( TF1* f, bool valid, double chisquare ):
    _function( f ),
    _valid( valid ),
    _chisquare( chisquare )
  {}

  constexpr bool valid() const
  { return _function && _valid && _chisquare > 0; }

  TF1* _function = nullptr;
  bool _valid = false;
  double _chisquare = -1;

};

//_______________________________________
inline bool operator < (const Result& lhs, const Result& rhs )
{

  if( !(lhs._function && lhs._valid && lhs._chisquare > 0 ) ) return false;
  else if( !(rhs._function && rhs._valid && rhs._chisquare > 0 ) ) return true;
  else return lhs._chisquare < rhs._chisquare;

}

template<class T>
bool IsSuccess( T fitResult )
{ return fitResult.Get() && fitResult->Status() == 0 && fitResult->CovMatrixStatus() == 3; }

//_______________________________________
Double_t FitFunction( Double_t* x, Double_t* par )
{ return par[0]*std::exp( -0.5*square( (*x-par[1])/par[2] ) ); }

//_______________________________________
Result Fit( TH1* h )
{
  if( !h->GetEntries() ) return {};

  const auto fname = Form( "%s_fit", h->GetName() );
  const auto xMin = h->GetXaxis()->GetXmin();
  const auto xMax = h->GetXaxis()->GetXmax();

//   auto FitFunction = []( Double_t* x, Double_t* par )
//   { return par[0]*std::exp( -0.5*square( (*x-par[1])/par[2] ) ); };

  auto f = new TF1( fname, FitFunction, xMin, xMax, 3 );
  f->SetNpx(500);
  f->SetParameter(0, h->GetMaximum() );
  f->SetParameter(1, h->GetMean() );
  f->SetParameter(2, h->GetRMS()/2 );

  f->SetLineColor(2);
  f->SetLineWidth(2);

  const auto fitResult = h->Fit( f, "0QS" );

  return { f, IsSuccess( fitResult ), f->GetChisquare()/f->GetNDF() };

}

//_______________________________________
Double_t FitFunction2( Double_t* x, Double_t* par )
{ 
  return 
    par[0]*std::exp( -0.5*square( (*x-par[1])/par[2] ) ) +
    par[3]*std::exp( -0.5*square( (*x-par[1])/par[4] ) );
}

//_______________________________________
Result Fit_double( TH1* h )
{
  if( !h->GetEntries() ) return {};

  const auto fname = Form( "%s_fit", h->GetName() );
  const auto xMin = h->GetXaxis()->GetXmin();
  const auto xMax = h->GetXaxis()->GetXmax();

  auto f = new TF1( fname, FitFunction2, xMin, xMax, 5 );
  f->SetNpx(500);
  f->SetParameter(0, h->GetMaximum() );
  f->SetParameter(1, h->GetMean() );
  f->SetParameter(2, h->GetRMS()/2 );
  f->SetParameter(3, h->GetMaximum()/10 );
  f->SetParameter(4, h->GetRMS()*2 );

  f->SetLineColor(2);
  f->SetLineWidth(2);

  const auto fitResult = h->Fit( f, "0S" );

  return { f, IsSuccess( fitResult ), f->GetChisquare()/f->GetNDF() };

}

//_______________________________________
Double_t FitFunction3( Double_t* x, Double_t* par )
{
  const Double_t t0 = (x[0]-par[1]+par[2]*std::sqrt(12)/2)/par[3];
  const Double_t t1 = (x[0]-par[1]-par[2]*std::sqrt(12)/2)/par[3];
  return par[0]*( std::erf(t0) + std::erf(-t1) );
}

//_______________________________________
Result Fit_box( TH1* h )
{

  if( !h->GetEntries() ) return {};

  const auto fname = Form( "%s_fit", h->GetName() );
  const auto xMin = h->GetXaxis()->GetXmin();
  const auto xMax = h->GetXaxis()->GetXmax();

  auto f = new TF1( fname, FitFunction3, xMin, xMax, 4 );
  f->SetNpx(500);
  f->SetParameter(0, h->GetMaximum() );
  f->SetParameter(1, h->GetMean() );
  f->SetParameter(2, h->GetRMS() );
  f->SetParameter(3, 0.5 );

  f->SetLineColor(2);
  f->SetLineWidth(2);

  TFitResultPtr fitResult = nullptr;
  for( int iTry = 0; iTry < 10; ++iTry )
  {

    fitResult =  h->Fit( f, "0SQ" );
    if( fitResult->Status() == 0 && fitResult->CovMatrixStatus() == 3 ) break;

  }

  /* we add a 1.5 factor to the chisquare to slightly favor gaussian fit over box fit */
  return { f, IsSuccess( fitResult ), 1.5*f->GetChisquare()/f->GetNDF() };

}
#endif
