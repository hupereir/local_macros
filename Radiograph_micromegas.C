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

namespace
{

  double get_sector_phi( int isec ) 
  { return isec*M_PI/6; }

  double normalize_phi( double phi ) 
  { return phi > M_PI ? phi - 2*M_PI:phi; }
  
}

//____________________________________________________________________________
void Radiograph_micromegas()
{

  set_style( false );

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // input files
  const TString tag = "_acts_full_notpc_nodistortion";
  const TString inputFile = Form( "DST/CONDOR_realistic_micromegas/dst_reco%s/dst_reco*_?.root", tag.Data() );
  
  const TString pdfFile = Form( "Figures/Radiograph_micromegas%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  if( true )
  {
    
    const TString var( "_clusters._truth_z:_clusters._truth_x" );
    
    // truth information
    for( int i = 0; i < 2; ++i )
    {
      
      const TCut layer_cut = Form( "_clusters._layer==%i", 55+i );
      
      // create canvas
      auto cv = new TCanvas( "cv", "cv",900, 900 );
           
      const TString hname = Form( "radiograph_truth_%i", 55+i );
      auto h = new TH2F( hname, hname, 500, -100, 100, 500, -105, 105 );
      Utils::TreeToHisto( tree, hname, var, layer_cut, false );
      
      h->GetXaxis()->SetTitle( "x (cm)" );
      h->GetYaxis()->SetTitle( "z (cm)" );
      h->SetTitle( hname );
      
      h->Draw();
      gPad->Update(); 
      
      Draw::VerticalLine( gPad, 0 )->Draw();
      
      pdfDocument.Add( cv );
      
    }
  }
  
  if( false )
  {
    
    const TString var( "_clusters._truth_z:_clusters._truth_phi" );
    
    // truth information
    for( int i = 0; i < 2; ++i )
    {
      
      const TCut layer_cut = Form( "_clusters._layer==%i", 55+i );
      
      // create canvas
      auto cv = new TCanvas( "cv", "cv", 900, 300 );
      
      cv->Divide( 3, 1 );
      
      const TString hname = Form( "radiograph_truth_%i", 55+i );
      auto h = new TH2F( hname, hname, 500, -3.14, 3.14, 500, -105, 105 );
      Utils::TreeToHisto( tree, hname, var, layer_cut, false );
      
      h->GetXaxis()->SetTitle( "#phi (rad)" );
      h->GetYaxis()->SetTitle( "z (cm)" );
      h->SetTitle( hname );
      
      cv->cd(1);
      h->Draw();
      gPad->Update(); 
      
      // draw sector boundaries
      for( int i = 0; i < 12; ++i )
      { 
        Draw::VerticalLine( gPad, normalize_phi(get_sector_phi(i)+M_PI/12 ))->Draw(); 
      }
      
      cv->cd(2);
      h->ProjectionX()->Draw();
      gPad->Update(); 
      
      // draw sector boundaries
      for( int i = 0; i < 12; ++i )
      { Draw::VerticalLine( gPad, normalize_phi(get_sector_phi(i)+M_PI/12 ))->Draw(); }
      
      cv->cd(3);
      h->ProjectionY()->Draw();
      
      pdfDocument.Add( cv );
      
    }
  }
  
  // create canvas
  if( false )
  {
    auto cv = new TCanvas( "cv_2", "cv_2", 900, 300 );
    cv->Divide( 2, 1 );

    {
      const TString var( "_clusters._phi" );
      const TCut layer_cut( "_clusters._layer==55");
  
      const TString hname( "radiograph_55" );
      auto h = new TH1F( hname, hname, 500, -3.14, 3.14 );
      Utils::TreeToHisto( tree, hname, var, layer_cut, false );
  
      h->GetXaxis()->SetTitle( "#phi (rad)" );
      h->SetTitle("");
      cv->cd(1);
      h->Draw();
    }
    

    if(false )
    {
      const TString var( "_clusters._z" );
      const TCut layer_cut( "_clusters._layer==56");
  
      const TString hname( "radiograph_56" );
      auto h = new TH1F( hname, hname, 500, -105, 105 );
      Utils::TreeToHisto( tree, hname, var, layer_cut, false );
  
      h->GetXaxis()->SetTitle( "z (cm)" );
      h->SetTitle("");
      cv->cd(2);
      h->Draw();
    }
 
    pdfDocument.Add( cv );
  
  }
    
}
