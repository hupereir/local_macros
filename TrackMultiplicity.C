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
TString TrackMultiplicity()
{
  
  set_style( false );
  
  const TString tag = "_Hijing_Micromegas_0_20fm";
  const TString inputFile = "DST/CONDOR_Hijing_Micromegas/dst_eval_sHijing_0_20fm-new/dst_eval_*.root";

  const TString pdfFile = Form( "Figures/TrackMultiplicity%s.pdf", tag.Data() );
  std::cout << "TrackMultiplicity - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );
  
  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  
  {
    // tracks per events
    auto cv( new TCanvas( "cv0", "cv0", 800, 800 ) );
    cv->SetRightMargin( 0.1 );
    const TString hname = "h";
    auto h = new TH1F( hname, hname, 100, 0, 1500 );
    const TString var = "DST#EVAL#SimEvaluator_hp::Container._events[0]._nparticles";
    const TCut cut = "DST#EVAL#SimEvaluator_hp::Container._events[0]._bimp < 15";
    Utils::TreeToHisto( tree, hname, var, cut, false );

    h->SetTitle("");
    h->GetXaxis()->SetTitle( "#tracks/event" );
    h->GetXaxis()->SetMaxDigits(3);
    h->Draw();
    
    cv->SetLogy( true );

    Draw::PutText( 0.5, 0.8, Form( "mean: %.0f #pm %.0f", h->GetMean(), h->GetMeanError() ) );
    
    pdfDocument.Add( cv );   
  }

  TProfile* p = nullptr; 
  {
    // tracks per events vs centrality
    auto cv( new TCanvas( "cv1", "cv1", 800, 800 ) );
    cv->SetLeftMargin( 0.18 );

    // track multiplicity vs centrality
    const TString pname = "profile";
    p = new TProfile( pname, pname, 20, 0, 100 );
    const TString var = 
      "DST#EVAL#SimEvaluator_hp::Container._events[0]._nparticles:"
      "get_centrality(DST#EVAL#SimEvaluator_hp::Container._events[0]._bimp)";
    
    const TCut cut = "DST#EVAL#SimEvaluator_hp::Container._events[0]._bimp < 15";

    Utils::TreeToHisto( tree, pname, var, cut, false );

    p->SetTitle("");
    p->GetXaxis()->SetTitle( "Centrality (%)" );
    p->GetYaxis()->SetTitle( "#tracks/event" );
    p->GetYaxis()->SetTitleOffset( 1.8 );
      
    p->SetMarkerStyle(20);
    p->Draw();
    
    pdfDocument.Add( cv );   
  }
  
  if( p )
  {
    // tracks per events vs centrality
    auto cv( new TCanvas( "cv2", "cv2", 800, 800 ) );
    cv->SetLeftMargin( 0.15 );

    // convert to histogram
    auto h = new TH1F( "h0", "h0", p->GetXaxis()->GetNbins(), p->GetXaxis()->GetXmin(), p->GetXaxis()->GetXmax() );
    h->GetXaxis()->SetTitle( p->GetXaxis()->GetTitle() );
    h->GetYaxis()->SetTitle( "#tracks (c>c_{0}) / #tracks_{tot}" );
    for( int i =0; i < p->GetXaxis()->GetNbins(); ++i )
    { 
      h->SetBinContent( i+1, p->GetBinContent( i+1 ) );
      h->SetBinError( i+1, p->GetBinError( i+1 ) );
    }
    
    // integrated tracks per event
    auto hi = Utils::Integrate( h, true, true );
    hi->SetTitle( "" );
    hi->SetMarkerStyle(20);
    hi->Draw();
    pdfDocument.Add( cv );   
  }

  return pdfFile;
      
}
