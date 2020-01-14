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
void InvMassUpsilon( TString tag = TString() )
{

  set_style( false );

  //if( tag.IsNull() ) tag = "_5k_upsilon_full";
  // if( tag.IsNull() ) tag = "_5k_upsilon_truth";
  if( tag.IsNull() ) tag = "_5k_upsilon_truth_noouter";
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );
  const TString pdfFile = Form( "Figures/InvMassUpsilon%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  const TString var( "_track_pairs[]._m" );
  const TCut charge_cut( "_track_pairs[]._charge==0" );
  const TCut mask_cut( "HasOuterTracker(_track_pairs[]._trk_mask[0]) && HasOuterTracker(_track_pairs[]._trk_mask[1])" );
  const TCut momentum_cut( "_track_pairs[]._trk_pt[0]>6 && _track_pairs[]._trk_pt[1]>6" );

  const TString hname( "invMassH" );
  auto h = new TH1F( hname, "", 100, 5, 12 );
  Utils::TreeToHisto( tree, hname, var, charge_cut&&mask_cut, false );
  h->SetTitle( "" );

  std::cout << "InvMassUpsilon - Entries: " << h->GetEntries() << std::endl;

  auto cv = new TCanvas( "cv", "cv", 800, 800 );
  h->SetMarkerStyle( 20 );
  h->GetXaxis()->SetTitle( "M_{e+e-} (GeV/#it{c}^{2})" );
  h->SetMaximum( 1.2*h->GetMaximum() );
  h->Draw( "E" );

  const auto binWidth = h->GetXaxis()->GetBinWidth(1);
  std::cout << "InvMassUpsilon - bin width: " << binWidth << std::endl;
  static constexpr float mass_upsilon_1s = 9.460 ;

  auto f = new TF1( "invMassF", UTILS::FitUtils::CrystallBall, 5, 12, 5 );
  f->SetNpx(1000);
  f->SetParameter(0, binWidth*h->GetEntries() );
  f->SetParameter(1, mass_upsilon_1s );
  f->SetParameter(2, 0.050 );
  f->SetParameter(3, 1 );
  f->SetParameter(4, 3 );

  f->SetParLimits( 4, 1.05, 10 );

  h->Fit( f, "0I" );
  f->SetLineColor( 2 );
  f->Draw( "same" );

  auto text = new TPaveText( 0.16, 0.8, 0.5, 0.9, "NDC" );
  text->SetFillColor(0);
  text->SetFillStyle(0);
  text->SetBorderSize(0);
  text->SetTextAlign(11);

  text->AddText( Form( "Counts: %.0f #pm %.0f", f->GetParameter(0)/binWidth, f->GetParError(0)/binWidth ) );
  text->AddText( Form( "Mass: %.3f #pm %.3f GeV/#it{c}^{2}", f->GetParameter(1), f->GetParError(1) ) );
  text->AddText( Form( "Width: %.1f #pm %.1f MeV/#it{c}^{2}", 1000*f->GetParameter(2), 1000*f->GetParError(1) ) );

  Draw::UpdatePaveSize( text );
  text->Draw();

  double integralError = 0;
  double integral = h->IntegralAndError( 1, h->GetNbinsX(), integralError );

  std::cout << "InvMassUpsilon -"
    << " bins: " << h->GetNbinsX()
    << " Integral: " << integral << "+/-" << integralError
    << std::endl;

  pdfDocument.Add( cv );
}
