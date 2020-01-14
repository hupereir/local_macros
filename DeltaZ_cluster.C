#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
TString DeltaZ_cluster( TString tag = TString() )
{

  set_style( false );
  constexpr std::array<Float_t, nDetectors> maxDetResidual = { 0.003, 1, 0.3, 0.3, 0.3, 0.5 };
  auto maxResidual = *std::max_element( maxDetResidual.cbegin(), maxDetResidual.cend() )/3;

  // pdf output
  if( tag.IsNull() ) tag = "_1k_full_notpc_nz500" ;
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );
  const TString pdfFile = Form( "Figures/DeltaZ_cluster%s.pdf", tag.Data() );

  std::cout << "DeltaZ - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaZ - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool doFit = false;
  const int firstBoxLayer = 55;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "_clusters._z - _clusters._truth_z" );
  const TString var2d = Form( "%s:_clusters._layer", var.Data() );
  // const TCut cluster_cut( "_clusters._size == 1" );
  const TCut cluster_cut;

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {
    const TString hname( Form( "deltaz_%i", idet ) );
    std::unique_ptr<TH2> h2d( new TH2F( hname, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 100, -maxDetResidual[idet], maxDetResidual[idet] ) );
    Utils::TreeToHisto( tree, hname, var2d, cluster_cut, false );

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    auto cv = new TCanvas( cvName, cvName, 800, 800 );
    Draw::DivideCanvas( cv, nLayers[idet], false );

    // loop over layers
    for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {

      int layerIndex = firstLayer[idet] + ilayer;
      const auto hname = Form( "h_%i", layerIndex );
      TH1* h = h2d->ProjectionY( hname, ilayer+1, ilayer+1 );
      h->SetTitle( hname );
      h->SetLineColor( 1 );
      h->SetMarkerColor( 1 );
      h->GetXaxis()->SetTitle( "#Deltaz_{cluster-truth} (cm)" );
      h->GetXaxis()->SetRangeUser( -maxDetResidual[idet], maxDetResidual[idet] );
      h->SetMaximum( h->GetMaximum()*1.2 );

      cv->cd( ilayer+1 );
      h->Draw();

      // fit
      if( h->GetEntries() )
      {
        if( doFit )
        {
          auto f = (firstBoxLayer>=0 && layerIndex>=firstBoxLayer) ? Fit_box( h ):Fit( h );
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

      // draw vertical line at zero
      gPad->Update();
      Draw::VerticalLine( gPad, 0 )->Draw();

    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv );

  }

  return pdfFile;

}
