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
#include <TStyle.h>

#include "LayerDefines.h"
#include "Fit.C"

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
Float_t delta_phi( Float_t phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
TString DeltaRPhi_cluster_charge( TString tag = TString() )
{

  set_style( false );

  // initial guess for max residuals
  std::array<Float_t, nDetectors> max_det_residual = { 0.01, 0.01, 1.2, 1.2, 1.2, 1.5};

  if( tag.IsNull() ) tag = "_5k_flat_full_nominal_new" ;
  const TString inputFile = Form( "DST/dst_reco%s.root", tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaRPhi_cluster%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/DeltaRPhi_cluster%s.root", tag.Data() );

  std::cout << "DeltaRPhi_cluster - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhi_cluster - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "_tracks._clusters._r*delta_phi(_tracks._clusters._phi - _tracks._clusters._truth_phi)" );
  const TString var2d = Form( "%s:_tracks._clusters._layer", var.Data() );
  const TCut cluster_cut;
  const TCut momentum_cut;

  static constexpr int ncuts = 2;
  const std::array<TCut, 2> charge_cut = 
  {{
    "_tracks._charge < 0", 
    "_tracks._charge > 0" 
  }}; 
      
  // labels
  static const std::array<TString,2> label = {{ "q<0", "q>0" }};
  static constexpr std::array<int,2> color = {{ 1, 2 }};
  
  // optimize max residual
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    const TString hname( Form( "deltarphi_%i_0", idet ) );
    const TCut layer_cut( Form( "_tracks._clusters._layer==%i", firstLayer[idet]+1 ) );

    for( int i=0; i<3; ++i )
    {
      std::unique_ptr<TH1> h1( new TH1F( hname, "", 500, -max_det_residual[idet], max_det_residual[idet] ) );
      Utils::TreeToHisto( tree, hname, var, momentum_cut&&layer_cut, false );
      max_det_residual[idet] = 5*h1->GetRMS() + std::abs(h1->GetMean() );

    }

  }

  // create TGraph to store resolution vs layer
  auto tg = new TGraphErrors();
  tg->SetName( "residuals" );

  // create TGraph to store resolution vs radius
  auto tgl = new TGraphErrors();
  tgl->SetName( "residuals_layers" );

  // save all histograms
  std::array<TH1*, nLayersTotal*ncuts> h_array;

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, (idet == 0 ) ? 400:800 ) );
    Draw::DivideCanvas( cv.get(), nLayers[idet], false );
    
    for( int icut = 0; icut < ncuts; ++icut )
    {
    
      const TString hname( Form( "deltarphi_%i_%i", idet, icut ) );
      std::unique_ptr<TH2> h2d( new TH2F( hname, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 100, -max_det_residual[idet], max_det_residual[idet] ) );
      Utils::TreeToHisto( tree, hname, var2d, cluster_cut&&momentum_cut&&charge_cut[icut], false );

      // loop over layers
      for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
      {
        
        int layerIndex = firstLayer[idet] + ilayer;
        const auto hname = Form( "h_%i_%i", layerIndex, icut );
        TH1* h = h2d->ProjectionY( hname, ilayer+1, ilayer+1 );
        h->SetTitle( hname );
        h->SetLineColor( color[icut] );
        h->SetMarkerColor( 1 );
        h->GetXaxis()->SetTitle( "r.#Delta#phi_{clus-truth} (cm)" );
        h->GetXaxis()->SetRangeUser( -max_det_residual[idet], max_det_residual[idet] );
        h->SetMaximum( h->GetMaximum()*1.2 );
        
        cv->cd( ilayer+1 );
        
        if( icut ) h->Draw( "same" );
        else h->Draw();
      
        // save in array
        h_array[ncuts*layerIndex+icut] = h;

        // draw vertical line at zero
        if( icut == 0 )
        {
          gPad->Update();
          Draw::VerticalLine( gPad, 0 )->Draw();
        }
        
      }
    
    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );

  }

  // TGraph
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtgl", "cvtgl", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, nLayersTotal ) );
    h->SetMinimum(0);
    h->SetMaximum(1.2*Utils::GetMaximum( tgl ));
    h->GetXaxis()->SetTitle( "layer id" );
    h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (cluster-truth) (#mum)" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();

    tgl->SetMarkerStyle(20);
    tgl->SetLineColor(1);
    tgl->SetMarkerColor(1);
    tgl->Draw("P");

    pdfDocument.Add( cv.get() );
  }

  // TGraph
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtg", "cvtg", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, 90 ) );
    h->SetMinimum(0);
    h->SetMaximum(1.2*Utils::GetMaximum( tg ));
    h->GetXaxis()->SetTitle( "r (cm)" );
    h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (cluster-truth) (#mum)" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();

    tg->SetMarkerStyle(20);
    tg->SetLineColor(1);
    tg->SetMarkerColor(1);
    tg->Draw("P");

    pdfDocument.Add( cv.get() );
  }

  // save everything in rootfiles
  std::unique_ptr<TFile> output( TFile::Open( rootFile, "RECREATE" ) );
  output->cd();
  for( auto&& h:h_array) { h->Write(); }
  tgl->Write();
  tg->Write();
  output->Close();

  return pdfFile;

}
