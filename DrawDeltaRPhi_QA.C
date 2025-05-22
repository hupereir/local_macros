#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

#include <TFile.h>
#include <TGraphErrors.h>
#include <TStyle.h>

#include "LayerDefines.h"
#include "Fit.C"

R__LOAD_LIBRARY(libRootUtilBase.so)

//__________________________________________
TString DrawDeltaRPhi_QA()
{

  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // configuration
  const bool do_fit = true;

  const TString file = "QA/qa_output.root";
  const TString pdfFile = "Figures/DeltaRPhi_QA.pdf";
  PdfDocument pdfDocument( pdfFile );

  std::unique_ptr<TFile> input( TFile::Open( file ) );

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    Draw::DivideCanvas( cv.get(), nLayers[idet], false );

    // loop over layers
    for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {

      int layerIndex = firstLayer[idet] + ilayer;
      const auto hname = Form( "h_QAG4SimulationMvtx_drphi_%i", layerIndex );
      TH1* h = static_cast<TH1*>( input->Get( hname ) );

      if( !h ) continue;

      cv->cd( ilayer+1 );
      h->Draw( "hist" );

      // fit
      if( h->GetEntries() )
      {
        if( do_fit )
        {
          const auto result = std::min( Fit( h ), Fit_box( h ) );
          auto f = result._function;
          f->Draw("same");
          auto h = f->GetHistogram();
          auto rms = h->GetRMS();
          auto error = f->GetParError(2);

          Draw::PutText( 0.2, 0.8, Form( "#sigma = %.3g #pm %.3g #mum", rms*1e4, error*1e4 ) );

        } else {

          auto rms = h->GetRMS();
          auto error = h->GetRMSError(2);
          Draw::PutText( 0.2, 0.8, Form( "#sigma = %.3g #pm %.3g #mum", rms*1e4, error*1e4 ) );

        }
      }
    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );

  }

  return pdfFile;

}
