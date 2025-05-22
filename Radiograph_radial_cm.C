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

float get_phi( float phi )
{
  while( phi < 0 ) phi += 2.*M_PI;
  while( phi >= 2.*M_PI ) phi -= 2*M_PI;
  return phi;
}

//____________________________________________________________________________
void Radiograph_radial_cm()
{

  set_style( false );

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

//   // input files
//   const TString tag = "_all_directlasers_simple_no_distortion-new";
//   const TString inputFile = Form( "DST/dst_reco%s.root", tag.Data() );

  // input files
  const TString tag = "_centralmembrane-nominal";
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  std::cout << "Radiograph - inputFile: " << inputFile << std::endl;
  
  // pdf output
  const TString pdfFile = Form( "Figures/Radiograph_radial_cm%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  if( true )
  {
    const TString var( "_cm_clusters._r:get_phi( _cm_clusters._phi)" );

    const TCut side_cuts[] = 
    {
      "_cm_clusters._z > 0",
      "_cm_clusters._z < 0"      
    };
    
    const TCut cuts[] = 
    {
      "",
      "_cm_clusters._nclusters==1",
      "_cm_clusters._nclusters==2"
    };
    
    const TString labels[] = 
    { 
      "z>0",
      "z>0, size 1",
      "z>0, size 2",
      "z<0",
      "z<0, size 1",
      "z<0, size 2"
    };
    
    static constexpr int nsidecuts = 2;
    static constexpr int ncuts = 3;

    for( int iside = 0; iside < nsidecuts; ++iside )
    {
      
      const auto cvname = Form( "cv_%i", iside );
      auto cv = new TCanvas( cvname, cvname, 1200, 400 );
      cv->Divide( ncuts, 1 );
    
      for( int icut = 0; icut<ncuts; ++icut )
      {
        
        std::cout << "Radiograph_cm - side: " << side_cuts[iside] << " cut: " << cuts[icut] << std::endl;
        cv->cd(icut+1);
        
        const auto hname = Form( "h2d_%i%i", iside,icut );
        auto h2d = new TH2F( hname, "", 1000, 0, 2*M_PI, 1000, 20, 78 );
        Utils::TreeToHisto( tree, hname, var, side_cuts[iside]&&cuts[icut], kFALSE );
        
        h2d->GetXaxis()->SetTitle( "#phi (rad)" );
        h2d->GetYaxis()->SetTitle( "r (cm)" );
        h2d->SetTitle( "" );
        h2d->SetMarkerColor(1);
        h2d->SetLineColor(1);
        h2d->SetMarkerStyle(20);
        h2d->SetMarkerSize(0.2);
        h2d->Draw();
        
        // write cuts and entries
        Draw::PutText( 0.15, 0.15, Form( "%s - entries: %.0f", labels[icut+ncuts*iside].Data(), h2d->GetEntries() ) );
      }
      
      pdfDocument.Add( cv );
      
    }
  }

  if( true )
  {
    
    /// distortion correction grid size along phi
    int m_phibins = 36;
    
    static constexpr float m_phiMin = 0;
    static constexpr float m_phiMax = 2.*M_PI;
    
    /// distortion correction grid size along r
    int m_rbins = 16;
    
    static constexpr float m_rMin = 20; // cm
    static constexpr float m_rMax = 78; // cm
    
    const TString var( "_cm_clusters._r:get_phi(_cm_clusters._phi)" );

    const TCut cuts[] = 
    {
      "_cm_clusters._z > 0",
      "_cm_clusters._z < 0"
    };
    
    const TString labels[] = 
    { 
      "z>0",
      "z<0"
    };
    
    static constexpr int ncuts = 2;
    
    auto cv = new TCanvas( "cv2", "", 1000, 500 );
    cv->Divide( ncuts, 1 );
    
    for( int icut = 0; icut<ncuts; ++icut )
    {

      std::cout << "Radiograph_radial_cm - cut: " << cuts[icut] << std::endl;
      
      cv->cd( icut+1 );
      
      const auto hname = Form( "h2d_%i", icut );
      auto h2d = new TH2F( hname, "", m_phibins, m_phiMin, m_phiMax, m_rbins, m_rMin, m_rMax );
      Utils::TreeToHisto( tree, hname, var, cuts[icut], kFALSE );

      h2d->GetXaxis()->SetTitle( "#phi (rad)" );
      h2d->GetYaxis()->SetTitle( "r (cm)" );
      h2d->SetTitle( "" );
      h2d->Draw("colz");

      // write cuts and entries
      Draw::PutText( 0.15, 0.15, Form( "%s - entries: %.0f", labels[icut].Data(), h2d->GetEntries() ) );

      // write cuts and entries
      gPad->SetRightMargin( 0.24 );
    
    }
    
    pdfDocument.Add( cv );

  }
  if( true )
  {
    // y vs z
    const TString var( "_cm_clusters._x:_cm_clusters._z" );
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
    const TString var( "_cm_clusters._z" );
    const TString hName = "radiograph";
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
