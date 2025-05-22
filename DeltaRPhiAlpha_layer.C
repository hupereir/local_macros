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

#include "LayerDefines.h"
#include "Fit.C"

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
float delta_phi( float phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
TString DeltaRPhiAlpha_layer( TString tag = TString() )
{

  set_style( false );

  // input files
  if( tag.IsNull() ) tag = "_genfit_truth_notpc_nodistortion";
  // if( tag.IsNull() ) tag = "_acts_truth_notpc_nodistortion";
  const TString inputFile = Form( "DST/CONDOR_realistic_micromegas/dst_reco%s/dst_reco_realistic_micromegas_1?.root", tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaRPhiAlpha_layer%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/DeltaRPhiAlpha_layer%s.root", tag.Data() );

  std::cout << "DeltaRPhiAlpha_layer - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhiAlpha_layer - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) {
    std::cout << "DeltaRPhi - invalid tree" << std::endl;
    return pdfFile;
  }

  // variable names
  const TString var( "_tracks._clusters._trk_r*delta_phi(_tracks._clusters._trk_phi - _tracks._clusters._phi)" );
  const TString var2d = Form( "%s:tan(_tracks._clusters._trk_alpha)", var.Data() );
  const TString var3d = Form( "%s:_tracks._clusters._layer", var2d.Data() );
  
  const TCut momentum_cut = ("_tracks._pt>0.5" );
  const TCut pattern_cut(
    "_tracks._truth_pt>0.5"
    "&&_tracks._nclusters_mvtx>2"
    "&&_tracks._nclusters_intt>=2"
    "&&_tracks._nclusters_micromegas>=2"
    );

  const TCut layer_cut("_tracks._clusters[]._layer>=7 && _tracks._clusters[]._layer<55" );
  
  // histogram
  const TString hname = "h3d";
  const double max_alpha = 0.5;
  const double max_residuals = 0.5;
  auto h3d = new TH3F( hname, hname, 48, 7, 55, 100, -max_alpha, max_alpha, 100, -max_residuals, max_residuals );
  Utils::TreeToHisto( tree, hname, var3d, momentum_cut&&pattern_cut&&layer_cut, false );

  const auto entries = h3d->GetEntries();
  std::cout << "DeltaRPhiAlpha_layer - entries: " << entries << std::endl;
   
  // loop over detectors
  for( int idet = 2; idet < 5; ++idet )
  {
  
    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    Draw::DivideCanvas( cv.get(), nLayers[idet], false );

    // loop over layers
    for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {
      
        int layerIndex = firstLayer[idet] + ilayer;
        int layerBin = layerIndex - firstLayer[2] + 1;

        h3d->GetXaxis()->SetRange( layerBin, layerBin );
        TH2* h2d = static_cast<TH2*>( h3d->Project3D( "zy" ) );

        const TString hname = Form( "h2d_%i", layerBin );
        h2d->SetName( hname );
        cv->cd(ilayer+1);

        h2d->SetTitle( Form( "layer %i", layerIndex ) );
        h2d->GetXaxis()->SetTitle( "tan(#alpha)" );
        h2d->GetYaxis()->SetTitle( "r#Delta#phi (track - cluster) (cm)" );
        h2d->Draw();
        
        
        auto p = h2d->ProfileX();
        
        // p->SetMarkerStyle(20);
        p->SetLineColor(2);
        p->SetMarkerColor(2);
        p->Draw("same");

        gPad->Update();
        auto line = Draw::HorizontalLine( gPad, 0 );
        line->SetLineColor(4);
        line->Draw();
        
        const auto entries = h2d->GetEntries();
        std::cout << "DeltaRPhiAlpha_layer - layer: " << layerIndex << " entries: " << entries << std::endl;
        
      }
   
      cv->Update();
      cv->cd(0);
      pdfDocument.Add( cv.get() );
  }     
  
  return pdfFile;
}
