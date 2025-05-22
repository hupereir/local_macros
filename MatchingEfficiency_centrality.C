#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>
#include <g4eval/TrackingEvaluator_hp.h>

R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
// centrality from impact parameter
float get_centrality( float b )
{
  static constexpr float bmax=14.75;
  static constexpr float centmax = 92.2;
  return centmax*::square(b/bmax);
}

//____________________________________________________________________________
TString MatchingEfficiency_centrality()
{
  set_style( false );

//   const TString tag = "_sHijing_0_20fm_50kHz_bkg_0_20fm-tony";
//   const TString inputFile = "DST/CONDOR_hijing_micromegas-tony/dst_eval/dst_eval*.root";

  const TString tag = "_sHijing_0_20fm_50kHz_bkg_0_20fm-new";
  const TString inputFile = "DST/CONDOR_hijing_micromegas-new/dst_eval/dst_eval*.root";

  const bool use_strict = true;

  const TString pdfFile = use_strict ?
    Form( "Figures/MatchingEfficiency_centrality%s_strict.pdf", tag.Data() ):
    Form( "Figures/MatchingEfficiency_centrality%s.pdf", tag.Data() );

  const TString rootFilename = use_strict ?
    Form( "Rootfiles/MatchingEfficiency_centrality%s_strict.root", tag.Data() ):
    Form( "Rootfiles/MatchingEfficiency_centrality%s.root", tag.Data() );

  std::cout << "MatchingEfficiency_centrality - inputFile: " << inputFile << std::endl;
  std::cout << "MatchingEfficiency_centrality - pdfFile: " << pdfFile << std::endl;
  std::cout << "MatchingEfficiency_centrality - rootFilename: " << rootFilename << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) {
    std::cout << "CheckLayers - invalid tree" << std::endl;
    return pdfFile;
  }

  // output rootfile
  RootFile rootfile( rootFilename );

  // define centrality bins
//   static constexpr int ncentbins = 3;
//   std::array<double,10> centmin = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90 };
//   std::array<double,10> centmax = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
//   constexpr std::array<int, 10> color = { kBlue, kCyan+2, kGreen+2, kOrange+1, kRed, kBlue, kCyan+2, kGreen+2, kOrange+1, kRed };
  static constexpr int ncentbins = 5;
  std::array<double,5> centmin = { 0, 20, 40, 60, 80 };
  std::array<double,5> centmax = { 20, 40, 60, 80, 100 };
  constexpr std::array<int, 5> color = { kBlue, kGreen+2, kRed, kCyan+2, kOrange+1 };

  // loop over layers
  for( int i = 0; i < 2; ++i )
  {
    const int layer = 55+i;
    std::cout << "MatchingEfficiency_centrality - layer: " << layer << std::endl;

    const TString var = "_tracks._pt:get_centrality(DST#EVAL#SimEvaluator_hp::Container._events[0]._bimp)";
    const TCut basecut =
      "_tracks._embed == 1"
      "&&_tracks._pt > 0.5"
      "&&_tracks._nclusters_tpc > 20"
      "&&_tracks._contributors > 80";

    const TCut refcut = basecut && Form( "TrackingEvaluator_hp::has_layer( _tracks._truth_mask, %i )", layer );
    const TCut foundcut = use_strict ?
      (refcut && Form( "TrackingEvaluator_hp::has_layer( _tracks._correct_mask_strict, %i )", layer )):
      (refcut && Form( "TrackingEvaluator_hp::has_layer( _tracks._correct_mask, %i )", layer ));

    // get 2d ref histogram
    const auto refname = Form( "href_%i", layer );
    auto href_2d = new TH2F( refname, "", 100, 0, 100, 20, 0, 20 );
    Utils::TreeToHisto( tree, refname, var, refcut, false );
    rootfile.Add( href_2d );

    // found histogram
    const auto foundname = Form( "hfound_%i", layer );
    auto hfound_2d = new TH2F( foundname, "", 100, 0, 100, 20, 0, 20 );
    Utils::TreeToHisto( tree, foundname, var, foundcut, false );
    rootfile.Add( hfound_2d );

    // define histograms
    std::array<TH1*,ncentbins> href_array = {{}};
    std::array<TH1*,ncentbins> hfound_array = {{}};
    std::array<TH1*,ncentbins> heff_array = {{}};

    // loop over centrality bins
    for( int icent = 0; icent < ncentbins; ++icent )
    {

      // find relevant centrality bins
      int centbin_min = href_2d->GetXaxis()->FindBin( centmin[icent] );
      int centbin_max = href_2d->GetXaxis()->FindBin( centmax[icent] )-1;

      std::cout
        << "MatchingEfficiency_centrality -"
        << " centrality: (" << centmin[icent] << ", " << centmax[icent] << ") -"
        << " bins: (" << centbin_min << ", " << centbin_max << ")" << std::endl;

      // reference histogram
      const auto refname = Form( "href_%i_%.0f_%.0f", layer, centmin[icent], centmax[icent] );
      auto href = href_2d->ProjectionY( refname, centbin_min, centbin_max );
      href->SetMarkerStyle( 20 );
      href->SetMarkerColor( color[icent] );
      href->SetLineColor( color[icent] );
      href_array[icent] = href;
      rootfile.Add(href);

      // found histogram
      const auto foundname = Form( "hfound_%i_%.0f_%.0f", layer, centmin[icent], centmax[icent] );
      auto hfound = hfound_2d->ProjectionY( foundname, centbin_min, centbin_max );
      hfound->SetMarkerStyle( 20 );
      hfound->SetMarkerColor( color[icent] );
      hfound->SetLineColor( color[icent] );
      hfound_array[icent] = hfound;
      rootfile.Add(hfound);

      // efficiency histogram
      const auto effname = Form( "heff_%i_%.0f_%.0f", layer, centmin[icent], centmax[icent] );
      auto heff = new TH1F( effname, "", 20, 0, 20 );
      Utils::DivideHistograms( hfound, href, heff );
      heff->GetYaxis()->SetTitle( Form( "Matching efficiency, layer %i", layer ) );
      heff->SetMaximum(1);
      heff->SetMinimum(0);
      heff->SetMarkerStyle( 20 );
      heff->SetMarkerColor( color[icent] );
      heff->SetLineColor( color[icent] );

      heff_array[icent] = heff;
      rootfile.Add(heff);

    }

    // plot
    std::cout << "MatchingEfficiency_centrality - plotting" << std::endl;
    auto cv( new TCanvas( "cv", "cv", 800, 800 ) );
    for( int icent = 0; icent < ncentbins; ++icent )
    {

      auto heff = heff_array[icent];
      if( !heff ) { std::cout << "MatchingEfficiency_centrality - invalid histogram" << std::endl; }
      heff->SetMarkerStyle( 20 );
      heff->SetMarkerColor( color[icent] );
      heff->SetLineColor( color[icent] );
      if( icent == 0 ) heff->Draw();
      else heff->Draw("same");
    }

    pdfDocument.Add(cv);

  }


  return pdfFile;
}
