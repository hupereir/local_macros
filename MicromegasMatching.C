#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

TString MicromegasMatching()
{
    
//   const TString inputfile = "Rootfiles/PHMicromegasTpcTrackMatching-oldgeom.root";
//   const TString outputfile = "Figures/PHMicromegasTpcTrackMatching-oldgeom.pdf";

  const TString inputfile = "Rootfiles/PHMicromegasTpcTrackMatching-newgeom.root";
  const TString outputfile = "Figures/PHMicromegasTpcTrackMatching-newgeom.pdf";

  auto tfile = TFile::Open( inputfile );
  
  PdfDocument pdfdocument( outputfile );
  
  TCanvas* cv = new TCanvas( "cv", "cv", 900, 500 );
  cv->Divide( 2, 1 );
  
  {
    auto h = static_cast<TH1*>( tfile->Get( "rphi_0" ) );
    h->SetTitle( "r#Delta#phi (TPC - MM), layer 55" );
    cv->cd(1);
    h->Draw();
  }
  
  {
    auto h = static_cast<TH1*>( tfile->Get( "z_1" ) );
    h->SetTitle( "#Deltaz (TPC - MM), layer 56" );
    cv->cd(2);
    h->Draw();
  }

  pdfdocument.Add( cv );
  
  return outputfile;
  
}
