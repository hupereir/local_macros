#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
// centrality from impact parameter
float get_centrality( float b )
{
  static constexpr float bmax=14.75;
  static constexpr float centmax = 92.2;

//   static constexpr float bmax=12;
//   static constexpr float centmax = 100;
  return centmax*::square(b/bmax);
}

//____________________________________________________________________________
TString HitMultiplicity()
{

  set_style( false );

//   const TString tag = "_Hijing_Micromegas_0_12fm";
//   const TString inputFile = "DST/CONDOR_hijing_micromegas/dst_eval_sHijing_0_12fm/dst_eval_*.root";

  const TString tag = "_Hijing_Micromegas_0_20fm-new";
  const TString inputFile = "DST/CONDOR_hijing_micromegas/dst_eval_sHijing_0_20fm-new/dst_eval_*.root";

  const TString pdfFile = Form( "Figures/HitMultiplicity%s.pdf", tag.Data() );
  std::cout << "HitMultiplicity - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // number of tiles
  static constexpr int ntiles = 4;

  {
    // hits per event per tile
    std::unique_ptr<TCanvas> cv( new TCanvas( "cv", "cv", 800, 400 ) );
    cv->Divide( 2, 1 );

    for( int i = 0; i < 2; ++i )
    {
      const TString hname = Form( "hhits_%i", i );
      auto h = new TH1F( hname, "", 120, 0, 120 );
      const auto var = Form( "DST#EVAL#TrackingEvaluator_hp::Container._events[]._nhits[%i]/%i", 55+i, ntiles );
      const TCut cut = "DST#EVAL#SimEvaluator_hp::Container._events[0]._bimp < 15";
      Utils::TreeToHisto( tree, hname, var, cut, false );
      h->SetTitle( "" );
      h->GetXaxis()->SetTitle( Form( "N_{hits, Micromegas} (layer %i)", 55+i ) );

      cv->cd(i+1);
      h->Draw();

      gPad->SetLogy(true );

      Draw::PutText( 0.6, 0.8, Form( "Mean = %.1f #pm %.1f", h->GetMean(), h->GetMeanError() ) );
      Draw::PutText( 0.6, 0.75, Form( "Occ = %.1f %%", 100*h->GetMean()/256 ) );

      std::cout << "HitMultiplicity - layer " << i+55 << " " << h->GetMean() << " hits/event/tile" << std::endl;
    }

    pdfDocument.Add( cv.get() );

  }

  {
    // impact parameter distribution
    std::unique_ptr<TCanvas> cv( new TCanvas( "cv2", "cv2", 800, 800 ) );
    auto h = new TH1F( "bimp", "", 100, 0, 20 );
    const TString var = "DST#EVAL#SimEvaluator_hp::Container._events[0]._bimp";
    Utils::TreeToHisto( tree, h->GetName(), var, TCut(), false );

    h->SetTitle("");
    h->GetXaxis()->SetTitle( "b_{imp} (fm)" );
    h->Draw();
    pdfDocument.Add( cv.get() );
  }

  {
    // hits per event per tile
    std::unique_ptr<TCanvas> cv( new TCanvas( "cv1", "cv1", 800, 400 ) );
    cv->Divide( 2, 1 );

    for( int i = 0; i < 2; ++i )
    {
      // get multiplicity vs centrality
      TString pname = Form( "profile_%i", i );
      auto p = new TProfile( pname, pname, 10, 0, 100 );
      auto var = Form( "(DST#EVAL#TrackingEvaluator_hp::Container._events[0]._nhits[%i]/%i) * 100./256:get_centrality(DST#EVAL#SimEvaluator_hp::Container._events[0]._bimp)", 55+i, ntiles );
      Utils::TreeToHisto( tree, pname, var, TCut(), false );

      p->SetTitle("");
      p->GetXaxis()->SetTitle( "Centrality (%)" );
      p->GetYaxis()->SetTitle( Form( "Occupancy (%%) (layer %i)", i+55 ) );

      cv->cd(i+1);
      p->SetMarkerStyle(20);
      p->Draw();

    }

    pdfDocument.Add( cv.get() );

  }

  return pdfFile;

}
