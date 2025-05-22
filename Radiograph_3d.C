#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
void Radiograph_3d()
{

  set_style( false );

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

//   // input files
//   const TString tag = "_directlasers-1";
//   const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  const TString tag = "_flat_acts_full-newgeom4";
  const TString inputFile = Form( "DST/CONDOR%s/dst*1?.root", tag.Data() );

  const TString pdfFile = Form( "Figures/Radiograph_3d%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // configuration
  const Bool_t doFit = true;

  FileManager fileManager( inputFile );

  auto tree = fileManager.GetChain( "T" );

  constexpr int ncut = 3;
  const std::array<TCut, ncut> cut= 
  {{
    "_clusters._layer == 54",
    "_clusters._layer == 55",
    "_clusters._layer == 56"
  }};
  
  const std::array<int, ncut> color = {{ 1, 2, 2 }};
  
//   constexpr int ncut = 1;
//   const std::array<TCut, ncut> cut= 
//   {{
//     TCut()
//   }};
//   
//   const std::array<int, ncut> color = {{ 1 }};

  // variable names
  const TString var( "_clusters._truth_y:_clusters._truth_x:_clusters._truth_z" );
  // const TString var( "_clusters._y:_clusters._x:_clusters._z" );
  
  // create canvas
  auto cv = new TCanvas( "cv", "cv", 800, 800 );
  auto cv2 = new TCanvas( "cv2", "cv2", 800, 800 );
  for( int i = 0; i<ncut; ++i )
  {      
    const TString hName = Form( "radiograph_%i", i );
    auto h3d = new TH3F( hName, "", 500, -105, 105, 500, -90, 90, 500, -90, 90 );
    Utils::TreeToHisto( tree, hName, var, cut[i], false );
    
    h3d->GetXaxis()->SetTitle( "z (cm)" );
    h3d->GetYaxis()->SetTitle( "x (cm)" );
    h3d->GetZaxis()->SetTitle( "y (cm)" );
    h3d->SetTitle( "" );
    h3d->SetMarkerColor( color[i] );
    h3d->SetLineColor( color[i] );
    
    cv->cd();
    if( i == 0 ) h3d->Draw();
    else h3d->Draw("same");
    
    // draw 2d projection
    cv2->cd();
    if( i == 0 ) h3d->Project3D( "zy" )->Draw();
    else h3d->Project3D( "zy" )->Draw("same");
        
  }
  pdfDocument.Add( cv );
  pdfDocument.Add( cv2 );

}
