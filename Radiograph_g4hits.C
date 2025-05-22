#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TCut.h>
#include <TEllipse.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libg4eval.so)

#include "LayerDefines.h"
#include "Fit.C"

template< class T>
T normalize_phi( T phi )
{
  while( phi >= 2*M_PI ) phi -= 2*M_PI;
  while( phi < 0 ) phi += 2*M_PI;
  return phi;
}

namespace
{
  double get_sector_phi( int isec ) { return isec*M_PI/6 + M_PI/12; }
}

//____________________________________________________________________________
void Radiograph_g4hits()
{

  set_style( false );

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // input files
  const TString tag = "_all_directlasers";
  const TString inputFile = Form( "DST/dst_reco%s.root", tag.Data() );

//   // input files
//   const TString tag = "_centralmembrane-nominal";
//   const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  // pdf output
  const TString pdfFile = Form( "Figures/Radiograph_g4hits%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  if( true )
  {
    
    const TString var( "_g4hits._y:_g4hits._x" );

    const TCut cuts[] = 
    {
      "_g4hits._z > 0",
      "_g4hits._z < 0"
    };
    
    const TString labels[] = 
    { 
      "z>0",
      "z<0"
    };
    
    static constexpr int ncuts = 2;
    
    auto cv = new TCanvas( "cv", "", 1000, 500 );
    cv->Divide( ncuts, 1 );
    
    for( int icut = 0; icut<ncuts; ++icut )
    {

      std::cout << "Radiograph_g4hits - cut: " << cuts[icut] << std::endl;
      
      cv->cd( icut+1 );
      
      const auto hname = Form( "h2d_%i", icut );
      auto h2d = new TH2F( hname, "", 1000, -90, 90, 1000, -90, 90 );
      Utils::TreeToHisto( tree, hname, var, cuts[icut], kFALSE );

      h2d->GetXaxis()->SetTitle( "x (cm)" );
      h2d->GetYaxis()->SetTitle( "y (cm)" );
      h2d->SetTitle( "" );
      h2d->SetMarkerColor(1);
      h2d->SetLineColor(1);
      h2d->SetMarkerStyle(20);
      h2d->SetMarkerSize(0.2);
      h2d->Draw();

      // write cuts and entries
      Draw::PutText( 0.15, 0.15, Form( "%s - entries: %.0f", labels[icut].Data(), h2d->GetEntries() ) );
      
      // draw ellipse
      {
        auto ellipse = new TEllipse( 0, 0, 30, 30 );
        ellipse->SetFillStyle(0);
        ellipse->SetLineColor(2);
        ellipse->Draw();
      }
      
      {
        auto ellipse = new TEllipse( 0, 0, 20, 20 );
        ellipse->SetFillStyle(0);
        ellipse->SetLineColor(4);
        ellipse->Draw();
      }
    }
    
    pdfDocument.Add( cv );

  }

  if( true )
  {
    // y vs z
    const TString var( "_g4hits._x:_g4hits._z" );
    const TString hName = "radiograph";
    auto h2d = new TH2F( hName, "", 1000, -105.5, 105.5, 1000, -90, 90 );
    Utils::TreeToHisto( tree, hName, var, TCut(), kFALSE );

    auto cv = new TCanvas( "", "", 800, 800 );
    h2d->GetXaxis()->SetTitle( "z (cm)" );
    h2d->GetYaxis()->SetTitle( "x (cm)" );
    h2d->SetTitle( "" );
    h2d->Draw( "" );
    pdfDocument.Add( cv );
  }

  if( true )
  {
    // z
    const TString var( "_g4hits._z" );
    const TString hName = "radiograph";
    // auto h2d = new TH2F( hName, "", 1000, -105.5, 105.5, 1000, -90, 90 );
    auto h = new TH1F( hName, "", 1000, -5, 10 );
    Utils::TreeToHisto( tree, hName, var, TCut(), kFALSE );

    auto cv = new TCanvas( "", "", 800, 800 );
    cv->SetLeftMargin( .2 );
    h->GetXaxis()->SetTitle( "z (cm)" );
    h->SetTitle( "" );
    h->Draw( "" );

    pdfDocument.Add( cv );
  }

}
