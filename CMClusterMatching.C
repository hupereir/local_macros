#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
void CMClusterMatching()
{

  set_style( false );

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  const TString inputFile = "PHTpcCentralMembraneMatcher.root";
  const TString pdfFile = "Figures/PHTpcCentralMembraneMatcher.pdf";
  PdfDocument pdfDocument( pdfFile );
  
  TFile input( inputFile );
  const TString hname[] = {
    "hdr1_single",
    "hdr2_single",
    "hdr3_single",
    "hdr1_double",
    "hdr2_double",
    "hdr3_double" 
  };
    
  auto cv = new TCanvas( "cv", "", 1200, 800 );
  cv->Divide( 3, 2 );
  for( int i = 0; i < 6; ++i )
  {
    
    auto h = static_cast<TH1*>( input.Get(hname[i]) );
    h->GetXaxis()->SetTitle( "#Delta r (cm)" );
    cv->cd( i+1 );
    h->Draw();
  }
  
  pdfDocument.Add( cv );
}    
    
  
    
