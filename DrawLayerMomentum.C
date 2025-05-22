#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TFile.h>
#include <TGraphErrors.h>
#include <TStyle.h>

#include "LayerDefines.h"

R__LOAD_LIBRARY(libRootUtilBase.so)

//__________________________________________
void DrawLayerMomentum()
{

  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  constexpr int nFiles = 2;
  const std::array<TString, nFiles> file =
  {{
    "Rootfiles/Momentum_truth_realistic_truth_notpc_noouter.root",
    "Rootfiles/Momentum_truth_realistic_full_notpc_noouter.root"
  }};

  const TString pdfFile = "Figures/momentum_truth_realistic_notpc_noouter.pdf";

  PdfDocument pdfDocument( pdfFile );

  constexpr std::array<int, nFiles> color = {{ kBlack, kRed }};
  const std::array<TString, nFiles> label = {{ "truth", "full" }};

  std::array<std::unique_ptr<TFile>, nFiles> tfiles;
  for( int i = 0; i < nFiles; ++i )
  { tfiles[i].reset( TFile::Open( file[i] ) ); }

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    Draw::DivideCanvas( cv.get(), nLayers[idet], false );

    // loop over files
    for( int i = 0; i < nFiles; ++i )
    {

      // loop over layers
      for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
      {

        int layerIndex = firstLayer[idet] + ilayer;
        const auto hname = Form( "h_%i", layerIndex );
        TH1* h = static_cast<TH1*>( tfiles[i]->Get( hname ) );

        h->SetLineColor( color[i] );
        // h->SetMarkerColor( color[i] );
        h->GetXaxis()->SetTitle( "#it{p}_{T} (GeV/#it{c})" );

        // normalize
        h->Scale( 1.0/h->GetEntries() );

        cv->cd( ilayer+1 );
        if( i == 0 ) h->Draw( "h" );
        else h->Draw( "h same" );

        gPad->Update();

      }

    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );

  }

}
