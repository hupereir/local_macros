#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>
#include <g4eval/TrackingEvaluator_hp.h>

R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

TString MatchingPurity()
{
  set_style( false );

//   const TString tag = "_sHijing_0_20fm_50kHz_bkg_0_20fm";
//   const TString inputFile = "DST/CONDOR_hijing_micromegas-tony/dst_eval-new/dst_eval*.root";

  const TString tag = "_sHijing_0_20fm_50kHz_bkg_0_20fm-new";
  const TString inputFile = "DST/CONDOR_hijing_micromegas-new/dst_eval/dst_eval*.root";

//   const TString tag = "_sHijing_0_20fm_50kHz_bkg_0_20fm";
//   const TString inputFile = "DST/CONDOR_hijing_micromegas-tony/dst_eval/dst_eval*.root";

//   const TString tag = "_flat_full_acts";
//   // const TString tag = "_flat_truth_acts";
//   const TString inputFile = Form( "DST/CONDOR%s/dst_eval*.root", tag.Data() );

  const bool use_strict = true;

  const TString pdfFile = use_strict ?
    Form( "Figures/MatchingPurity%s_strict.pdf", tag.Data() ):
    Form( "Figures/MatchingPurity%s.pdf", tag.Data() );

  std::cout << "MatchingPurity - inputFile: " << inputFile << std::endl;
  std::cout << "MatchingPurity - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) {
    std::cout << "MatchingPurity - invalid tree" << std::endl;
    return pdfFile;
  }

  std::array<TH1*,2> href_array = {{}};
  std::array<TH1*,2> hfound_array = {{}};

  {

    // momentum distribution
    auto cv( new TCanvas( "cv", "cv", 1200, 600 ) );
    cv->Divide( 2, 1 );

    // loop over layers
    for( int i = 0; i < 2; ++i )
    {
      const int layer = 55+i;
      const TString var = "_tracks._pt";

      const TCut basecut =
        "_tracks._embed == 1"
        "&&_tracks._pt > 0.5"
        "&&_tracks._nclusters_tpc > 20";

//       const TCut basecut =
//         "_tracks._pt > 0.5"
//         "&&_tracks._nclusters_tpc > 20";

      const TCut refcut = basecut && Form( "TrackingEvaluator_hp::has_layer( _tracks._mask, %i )", layer );
      const TCut foundcut = use_strict ?
        refcut && Form( "TrackingEvaluator_hp::has_layer( _tracks._correct_mask_strict, %i )", layer ):
        refcut && Form( "TrackingEvaluator_hp::has_layer( _tracks._correct_mask, %i )", layer );

      // reference histogram
      const TString hname_ref = Form( "href_%i", i );
      auto href = new TH1F( hname_ref, "", 40, 0, 20 );
      Utils::TreeToHisto( tree, hname_ref, var, refcut, false );

      href->SetTitle( "" );
      href->SetLineColor(1);
      href->GetXaxis()->SetTitle( "#it{p}_{T} (GeV/#it{c})" );

      // found histogram
      const TString hname_found = Form( "hfound_%i", i );
      auto hfound = new TH1F( hname_found, "", 40, 0, 20 );
      Utils::TreeToHisto( tree, hname_found, var, foundcut, false );

      hfound->SetTitle( "" );
      hfound->SetLineColor(2);
      hfound->GetXaxis()->SetTitle( "#it{p}_{T} (GeV/#it{c})" );

      cv->cd(i+1);
      href->Draw();
      hfound->Draw("same");

      // store
      href_array[i] = href;
      hfound_array[i] = hfound;

    }

    pdfDocument.Add( cv );

  }

  {
    // matching efficiency vs momentum
    auto cv( new TCanvas( "cv", "cv", 1200, 600 ) );
    cv->Divide( 2, 1 );

    // loop over layers
    for( int i = 0; i < 2; ++i )
    {

      const int layer = 55+i;
      const TString hname_eff = Form( "heff_%i", i );
      auto heff = static_cast<TH1*>( hfound_array[i]->Clone( hname_eff ) );
      heff->Reset();

      Utils::DivideHistograms( hfound_array[i], href_array[i], heff );

      heff->GetYaxis()->SetTitle( Form( "Matching purity, layer %i", layer ) );

      cv->cd(i+1);
      gPad->SetLeftMargin( 0.17 );
      heff->GetYaxis()->SetTitleOffset( 1.7 );
      heff->SetMaximum(1);
      heff->SetMinimum(0);
      heff->SetMarkerStyle( 20 );
      heff->SetLineColor(1);
      heff->Draw();

    }

    pdfDocument.Add( cv );
  }

  return pdfFile;
}
