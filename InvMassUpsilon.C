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

//____________________________________________________________________________
TString InvMassUpsilon()
{

  static constexpr double mass_min = 6;
  static constexpr double mass_max = 12;
  static constexpr int nbins = 50*(mass_max - mass_min );

  set_style( false );
  const TString tag = "_upsilon_pythia8";
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s_?.root", tag.Data(), tag.Data() );

  static constexpr bool use_double_crystalball = true;
  const TString pdfFile = use_double_crystalball ?
    Form( "Figures/InvMassUpsilon_cb2%s.pdf", tag.Data() ):
    Form( "Figures/InvMassUpsilon_cb1%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  const TString rootfilename = Form( "Rootfiles/InvMassUpsilon%s.root", tag.Data() );
  RootFile rootfile( rootfilename );

  std::cout << "InvMassUpsilon - inputFile: " << inputFile << std::endl;
  std::cout << "InvMassUpsilon - pdfFile: " << pdfFile << std::endl;

  // Utils::max_entries = 1000;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return {};

  std::cout << "InvMassUpsilon - entries: " << tree->GetEntries() << std::endl;

  if( false )
  {
    // generator level invariant mass
    const TString var( "_track_pairs[]._truth_m" );

    // const TCut pattern_cut;
    const TCut mvtx_cut(
      "_track_pairs[]._tracks[0]._nclusters_mvtx>2"
      "&&_track_pairs[]._tracks[1]._nclusters_mvtx>2"
      );

    const TCut tpc_cut(
      "_track_pairs[]._tracks[0]._nclusters_tpc>30"
      "&&_track_pairs[]._tracks[1]._nclusters_tpc>30" );

    const TCut pattern_cut = mvtx_cut&&tpc_cut;

    const TString hname( "invMassH_truth" );
    auto h = new TH1F( hname, "", nbins, mass_min, mass_max );
    Utils::TreeToHisto( tree, hname, var, pattern_cut, false );
    h->SetTitle( "" );

    // draw
    auto cv = new TCanvas( "cv", "cv", 800, 800 );
    cv->SetLeftMargin( 0.15 );
    h->SetMarkerStyle( 20 );
    h->GetXaxis()->SetTitle( "M_{e+e-} (GeV/#it{c}^{2})" );
    h->SetMaximum( 1.2*h->GetMaximum() );
    // h->Draw( "E" );
    h->Draw( "" );

    gPad->SetLogy( true );

    pdfDocument.Add( cv );
    rootfile.Add( h );

  }

  {
    // reconstructed invariant mass
    const TString var( "_track_pairs[]._m" );

    // const TCut pattern_cut;
    const TCut mvtx_cut(
      "_track_pairs[]._tracks[0]._nclusters_mvtx>2"
      "&&_track_pairs[]._tracks[1]._nclusters_mvtx>2"
      );

    const TCut tpc_cut(
      "_track_pairs[]._tracks[0]._nclusters_tpc>30"
      "&&_track_pairs[]._tracks[1]._nclusters_tpc>30" );

    const TCut pattern_cut = mvtx_cut &&tpc_cut;


    const TString hname( "invMassH" );
    auto h = new TH1F( hname, "", nbins, mass_min, mass_max );
    Utils::TreeToHisto( tree, hname, var, pattern_cut, false );
    h->SetTitle( "" );

    std::cout << "InvMassUpsilon - Entries: " << h->GetEntries() << std::endl;

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

    std::cout << "InvMassUpsilon -"
      << " tag: " << tag
      << " bins: " << h->GetNbinsX()
      << " width: " << 1000*f->GetParameter(2) << "+/-" << 1000*f->GetParError(1) << "MeV"
      << " Entries: " << h->GetEntries()
      << " Integral: " << integral << "+/-" << integralError
      << std::endl;

    gPad->SetLogy( true );
    pdfDocument.Add( cv );
    rootfile.Add( h );
  }

  return pdfFile;
}
