#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>
#include <g4eval/TrackingEvaluator_hp.h>

R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"


//____________________________________________________________________________
TString DrawMatchingEfficiency_centrality()
{
  set_style( false );

  static constexpr int ncentbins = 5;
  std::array<double,5> centmin = { 0, 20, 40, 60, 80 };
  std::array<double,5> centmax = { 20, 40, 60, 80, 100 };
  constexpr std::array<int, 5> color = { kBlue, kGreen+2, kRed, kCyan+2, kOrange+1 };

//   static constexpr int ncentbins = 10;
//   std::array<double,10> centmin = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90 };
//   std::array<double,10> centmax = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
//   constexpr std::array<int, 10> color = { kBlue, kCyan+2, kGreen+2, kOrange+1, kRed, kBlue, kCyan+2, kGreen+2, kOrange+1, kRed };

  const TString tag = "_sHijing_0_20fm_50kHz_bkg_0_20fm-new";
  const TString inputFile = Form( "Rootfiles/MatchingEfficiency_centrality%s.root", tag.Data() );
  const TString pdfFile = Form( "Figures/MatchingEfficiency_centrality%s.pdf", tag.Data() );

  auto tfile = TFile::Open( inputFile );

  PdfDocument pdfDocument( pdfFile );

  // loop over layers
  for( int i = 0; i < 2; ++i )
  {
    const int layer = 55+i;
    std::cout << "DrawMatchingEfficiency_centrality - layer: " << layer << std::endl;

    auto cv( new TCanvas( "cv", "cv", 800, 800 ) );
    auto legend = new TLegend( 0.15, 0.2, 0.6, 0.45, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);

    #if false
    for( int icent = 0; icent < ncentbins; ++icent )
    {

      const auto effname = Form( "heff_%i_%.0f_%.0f", layer, centmin[icent], centmax[icent] );
      const auto heff = static_cast<TH1*>( tfile->Get(effname) );
      heff->SetMarkerStyle( 20 );
      heff->SetMarkerColor( color[icent] );
      heff->SetLineColor( color[icent] );
      if( icent == 0 ) heff->Draw();
      else heff->Draw("same");

      legend->AddEntry( heff, Form( "%.0f<cent<%.0f",  centmin[icent], centmax[icent] ) );
    }

    #else
    const auto refname = Form( "href_%i", layer );
    const auto href_2d = static_cast<TH2*>( tfile->Get(refname) );

    const auto foundname = Form( "hfound_%i", layer );
    const auto hfound_2d = static_cast<TH2*>( tfile->Get(foundname) );
    for( int icent = 0; icent < ncentbins; ++icent )
    {
      // find relevant centrality bins
      int centbin_min = href_2d->GetXaxis()->FindBin( centmin[icent] );
      int centbin_max = href_2d->GetXaxis()->FindBin( centmax[icent] )-1;

      std::cout
        << "DrawMatchingEfficiency_centrality -"
        << " centrality: (" << centmin[icent] << ", " << centmax[icent] << ") -"
        << " bins: (" << centbin_min << ", " << centbin_max << ")" << std::endl;

      // reference histogram
      const auto refname = Form( "href_%i_%.0f_%.0f", layer, centmin[icent], centmax[icent] );
      auto href = href_2d->ProjectionY( refname, centbin_min, centbin_max );

      // found histogram
      const auto foundname = Form( "hfound_%i_%.0f_%.0f", layer, centmin[icent], centmax[icent] );
      auto hfound = hfound_2d->ProjectionY( foundname, centbin_min, centbin_max );

      // efficiency histogram
      const auto effname = Form( "heff_%i_%.0f_%.0f", layer, centmin[icent], centmax[icent] );
      auto heff = new TH1F( effname, "", 20, 0, 20 );
      Utils::DivideHistograms( hfound, href, heff );
      heff->GetYaxis()->SetTitle( Form( "Matching efficiency, layer %i", layer ) );
      heff->GetXaxis()->SetTitle( "#it{p}_{T} (GeV/#it{c})" );
      heff->SetMaximum(1);
      heff->SetMinimum(0);
      heff->SetMarkerStyle( 20 );
      heff->SetMarkerColor( color[icent] );
      heff->SetLineColor( color[icent] );

      if( icent == 0 ) heff->Draw();
      else heff->Draw("same");

      legend->AddEntry( heff, Form( "%.0f-%.0f%%",  centmin[icent], centmax[icent] ) );

    }

    #endif

    legend->Draw();

    pdfDocument.Add(cv);

  }


  return pdfFile;
}
