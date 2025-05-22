#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/FitUtils.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TH1.h>
#include <TPaveText.h>
#include <TStyle.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

namespace
{
  bool exclude_signal = true;
  const double mass_signal_min = 7;
  const double mass_signal_max = 11;
}

//____________________________________________________________________________
double fit_function_background(  double* x, double* par )
{
  const double mass = x[0];
  if( exclude_signal && (mass > mass_signal_min && mass < mass_signal_max) )
  {
    TF1::RejectPoint();
    return 0;
  }

  const double amplitude = par[0];
  const double mean = par[1];
  const double sigma = par[2];
  const double slope = par[3];
  return amplitude*UTILS::FitUtils::VWG( mass, mean, sigma, slope );
}

//____________________________________________________________________________
double fit_function_all(  double* x, double* par )
{
  // first 7 parameters are for CB2
  const double signal = UTILS::FitUtils::CrystalBall2( x, par );

  const double mass = x[0];
  const double amplitude = par[7];
  const double mean = par[8];
  const double sigma = par[9];
  const double slope = par[10];
  const double background = amplitude*UTILS::FitUtils::VWG( mass, mean, sigma, slope );

  return signal + background;

}

//____________________________________________________________________________
void FitInvMassUpsilon()
{
  set_style( false );

  const TString tag = "_upsilon_pythia8_acts_full_no_distortion";
  const TString inputFile = Form( "Rootfiles/InvMassUpsilon%s.root", tag.Data() );
  FileManager fileManager( inputFile );

  const TString pdfFile = Form( "Figures/InvMassUpsilon%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  auto h = fileManager.GetHistogram( "invMassH" );
  assert(h);

  const double mass_min = h->GetXaxis()->GetXmin();
  const double mass_max = h->GetXaxis()->GetXmax();

  auto cv = new TCanvas( "cv", "cv", 800, 800 );
  cv->SetLeftMargin( 0.15 );
  h->SetMarkerStyle( 20 );
  h->GetXaxis()->SetTitle( "M_{e+e-} (GeV/#it{c}^{2})" );
  h->SetMaximum( 1.2*h->GetMaximum() );
  // h->SetMaximum( 3000 );
  h->Draw( "E" );

  const auto binWidth = h->GetXaxis()->GetBinWidth(1);
  std::cout << "InvMassUpsilon - bin width: " << binWidth << std::endl;
  static constexpr float mass_upsilon_1s = 9.460 ;

  auto f_bg = new TF1( "fit_function_background", fit_function_background, mass_min, mass_max, 4 );

  f_bg->SetNpx( 250*(mass_max - mass_min) );
  f_bg->SetParameter(0, binWidth*h->GetEntries()*1.2 );
  f_bg->SetParameter(1, 0 );
  f_bg->SetParameter(2, 1 );
  f_bg->SetParameter(3, 0.1 );

  // do the fit
  h->Fit( f_bg, "0I" );
  f_bg->SetLineColor(2);
  f_bg->Draw("same");

  if( true )
  {
    // print integrals:
    exclude_signal = false;
    const double total = binWidth*h->GetEntries();
    const double background = f_bg->Integral( mass_min, mass_max );
    std::cout << "FitInvMassUpsilon - total: " << total << std::endl;
    std::cout << "FitInvMassUpsilon - background: " << background << std::endl;

    // full signal
    auto f =  new TF1( "fit_function_all", fit_function_all, mass_min, mass_max, 11 );

    f->SetNpx( 250*(mass_max - mass_min) );
    f->SetParNames( "A_{#Upsilon(1S)}", "M", "#sigma", "#alpha_{1}", "n_{1}", "#alpha_{2}", "n_{2}" );
    f->SetParameter(0, (total - background)*1.2 );
    f->SetParameter(1, mass_upsilon_1s );
    f->SetParameter(2, 0.2 );
    f->SetParLimits( 2, 0.01, 0.8 );

    // low mass tail
    f->SetParameter(3, 0.7);
    f->SetParameter(4, 1.1);
    f->SetParLimits(3, 0.120, 10);
    f->SetParLimits(4, 1.05, 10);

    // high mass tail
    f->SetParameter(5, 2 );
    f->SetParameter(6, 1.1 );
    f->SetParLimits(5, 0.1, 10);
    f->SetParLimits(6, 1.05, 10);

    // background
    for( int i = 0; i < 4; ++i )
    // { f->SetParameter( i+7, f_bg->GetParameter(i) ); };
    { f->FixParameter( i+7, f_bg->GetParameter(i) ); };

    h->Fit( f, "0I" );

    f->SetLineColor( 2 );
    f->Draw( "same" );
  }

  gPad->SetLogy( true );
  pdfDocument.Add( cv );

}
