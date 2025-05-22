#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/FitUtils.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TH1.h>
#include <TPaveText.h>
#include <TStyle.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

//____________________________________________________________________________
bool HasOuterTracker( int64_t mask )
{ return (mask&(int64_t(1)<<55)) && (mask&(int64_t(1)<<56)); }

//____________________________________________________________________________
void InvMassUpsilon_QA( TString tag = TString() )
{

  set_style( false );

  if( tag.IsNull() ) tag = "_upsilon_acts_full_micromegas_nominal";
  const TString inputFile = Form( "DST/CONDOR%s/qa_*.root", tag.Data() );

  static constexpr bool use_double_crystalball = true;

  const TString pdfFile = use_double_crystalball ?
    Form( "Figures/InvMassUpsilon_QA_cb2%s.pdf", tag.Data() ):
    Form( "Figures/InvMassUpsilon_QA_cb1%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto h = fileManager.GetHistogram( "h_QAG4SimulationUpsilon_nReco_Pair_InvMassReco" );
  std::cout << "InvMassUpsilon_QA - Entries: " << h->GetEntries() << std::endl;

  auto cv = new TCanvas( "cv", "cv", 800, 800 );
  cv->SetLeftMargin( 0.15 );
  h->SetMarkerStyle( 20 );
  h->GetXaxis()->SetTitle( "M_{e+e-} (GeV/#it{c}^{2})" );
  h->GetXaxis()->SetRangeUser( 7, 11 );
  h->SetMaximum( 1.2*h->GetMaximum() );
  h->Draw( "E" );

  const auto binWidth = h->GetXaxis()->GetBinWidth(1);
  std::cout << "InvMassUpsilon_QA - bin width: " << binWidth << std::endl;
  static constexpr float mass_upsilon_1s = 9.460 ;

  auto f = use_double_crystalball ?
    new TF1( "invMassF", UTILS::FitUtils::CrystalBall2, mass_min, mass_max, 7 ):
    new TF1( "invMassF", UTILS::FitUtils::CrystalBall, mass_min, mass_max, 5 );

  f->SetNpx( 250*(mass_max - mass_min) );
  f->SetParNames( "A_{#Upsilon(1S)}", "M", "#sigma", "#alpha_{1}", "n_{1}", "#alpha_{2}", "n_{2}" );
  f->SetParameter(0, binWidth*h->GetEntries()*1.2 );
  f->SetParameter(1, mass_upsilon_1s );
  f->SetParameter(2, 0.08 );
  f->SetParLimits( 2, 0.01, 0.8 );

  // low mass tail
  f->SetParameter(3, 1);
  f->SetParameter(4, 3);
  f->SetParLimits(3, 0.120, 10);
  f->SetParLimits(4, 1.05, 10);

  // high mass tail
  if( use_double_crystalball )
  {
    f->SetParameter(5, 1 );
    f->SetParameter(6, 5 );
    f->SetParLimits(5, 0.1, 10);
    f->SetParLimits(6, 1.05, 10);
  }

  // do the fit
  h->Fit( f, "0I" );
  f->SetLineColor( 2 );
  f->Draw( "same" );

  auto text = new TPaveText( 0.16, 0.8, 0.5, 0.9, "NDC" );
  // auto text = new TPaveText( 0.16, 0.65, 0.5, 0.9, "NDC" );
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(11);

  // text->AddText( Form( "Counts: %.0f #pm %.0f", f->GetParameter(0)/binWidth, f->GetParError(0)/binWidth ) );
  text->AddText( Form( "Mass: %.3f #pm %.3f GeV/#it{c}^{2}", f->GetParameter(1), f->GetParError(1) ) );
  text->AddText( Form( "Width: %.1f #pm %.1f MeV/#it{c}^{2}", 1000*f->GetParameter(2), 1000*f->GetParError(1) ) );

  if( false && use_double_crystalball )
  {
    const auto fractions = UTILS::FitUtils::CrystalBall2Fractions( f->GetParameter(3), f->GetParameter(4), f->GetParameter(5), f->GetParameter(6) );
    std::cout << "Fractions:" << std::endl;
    std::cout << Form( "  core: %.0f %%", fractions[1]*100 ) << std::endl;
    std::cout << Form( "  left tail: %.0f %%", fractions[0]*100 ) << std::endl;
    std::cout << Form( "  right tail: %.0f %%", fractions[2]*100 ) << std::endl;

    text->AddText( "Fractions:" );
    text->AddText( Form( "  core: %.0f %%", fractions[1]*100 ) );
    text->AddText( Form( "  left tail: %.0f %%", fractions[0]*100 ) );
    text->AddText( Form( "  right tail: %.0f %%", fractions[2]*100 ) );

    Draw::VerticalLine( cv, f->GetParameter(1)-f->GetParameter(2)*f->GetParameter(3) )->Draw();
    Draw::VerticalLine( cv, f->GetParameter(1)+f->GetParameter(2)*f->GetParameter(5) )->Draw();

  }

  text->Draw();
  Draw::UpdatePaveSize( text );

  double integralError = 0;
  double integral = h->IntegralAndError( 1, h->GetNbinsX(), integralError );

  std::cout << "InvMassUpsilon_QA -"
    << " tag: " << tag
    << " bins: " << h->GetNbinsX()
    << " width: " << 1000*f->GetParameter(2) << "+/-" << 1000*f->GetParError(1) << "MeV"
    << " Entries: " << h->GetEntries()
    << " Integral: " << integral << "+/-" << integralError
    << std::endl;

  pdfDocument.Add( cv );
}
