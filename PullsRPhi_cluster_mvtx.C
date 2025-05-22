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

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
Float_t delta_phi( Float_t phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
TString PullsRPhi_cluster_mvtx( TString tag = TString() )
{

  set_style( false );

  // initial guess for max residuals
  std::array<Float_t, nDetectors> max_det_residual = {{ 5, 5, 5, 5, 5, 5}};

  // input files
//   if( tag.IsNull() ) tag = "_realistic_full_nominal";
//   const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*_9*.root", tag.Data(), tag.Data() );

  if( tag.IsNull() ) tag = "_2k_realistic_full_zhengyun_fix3" ;
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  #if true
  const TString pdfFile = Form( "Figures/PullsRPhi_cluster%s_mvtx_0.pdf", tag.Data() );
  constexpr int n_cluster_cut = 3;
  std::array<TCut, n_cluster_cut> cluster_cuts =
  {{
    "_clusters._phi_size == 1",
    "_clusters._phi_size == 2",
    "_clusters._phi_size == 3"
  }};

  std::array<TString, n_cluster_cut> cut_label =
  {{
    "csize_{#phi} = 1",
    "csize_{#phi} = 2",
    "csize_{#phi} = 3"
  }};

  #else

  const TString pdfFile = Form( "Figures/PullsRPhi_cluster%s_mvtx_1.pdf", tag.Data() );
  constexpr int n_cluster_cut = 4;
  std::array<TCut, n_cluster_cut> cluster_cuts =
  {{
    "_clusters._phi_size == 2 && _clusters._z_size == 2 ",
    "_clusters._phi_size == 2 && _clusters._z_size == 3 ",
    "_clusters._phi_size == 3 && _clusters._z_size == 2 ",
    "_clusters._phi_size == 3 && _clusters._z_size == 3 "
  }};

  std::array<TString, n_cluster_cut> cut_label =
  {{
    "csize_{#phi} = 2, csize_{z} = 2",
    "csize_{#phi} = 2, csize_{z} = 3",
    "csize_{#phi} = 2, csize_{z} = 2",
    "csize_{#phi} = 3, csize_{z} = 3"
  }};

  #endif

  std::cout << "PullsRPhi_cluster - inputFile: " << inputFile << std::endl;
  std::cout << "PullsRPhi_cluster - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool do_fit = true;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // scale factor
  static const float pitch = 0.002688;
  static const float scale = pitch/std::sqrt(12.);

  // variable names
  const TString var = Form( "_clusters._r*delta_phi(_clusters._phi - _clusters._truth_phi)/%f", scale );
  // const TString var = "delta_phi(_clusters._phi - _clusters._truth_phi)/_clusters._phi_error";
  const TString var2d = Form( "%s:_clusters._layer", var.Data() );

  // optimize max residual
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    if( idet != 0 ) continue;

    const TString hname( Form( "PullsRPhi_%i_0", idet ) );
    const TCut layer_cut( Form( "_clusters._layer==%i", firstLayer[idet]+1 ) );

    for( int i=0; i<3; ++i )
    {
      std::unique_ptr<TH1> h1( new TH1F( hname, "", 500, -max_det_residual[idet], max_det_residual[idet] ) );
      Utils::TreeToHisto( tree, hname, var, layer_cut, false );
      max_det_residual[idet] = 3*h1->GetRMS() + std::abs(h1->GetMean() );
    }

  }

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {

    if( idet != 0 ) continue;

    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    cv->Divide( nLayers[idet], n_cluster_cut );

    // loop over cluster cuts
    for( int icut = 0; icut < n_cluster_cut; ++icut )
    {

      const TString hname( Form( "PullsRPhi_%i_%i", idet, icut ) );
      std::unique_ptr<TH2> h2d( new TH2F( hname, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 100, -max_det_residual[idet], max_det_residual[idet] ) );
      Utils::TreeToHisto( tree, hname, var2d, cluster_cuts[icut], false );

      // loop over layers
      for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
      {

        int layerIndex = firstLayer[idet] + ilayer;
        const auto hname = Form( "h_%i_%i", layerIndex, icut );
        TH1* h = h2d->ProjectionY( hname, ilayer+1, ilayer+1 );
        h->SetTitle( hname );
        h->SetLineColor( 1 );
        h->SetMarkerColor( 1 );
        h->GetXaxis()->SetTitle( "#Delta#phi_{clus-truth}/#sigma_{#phi}" );
        // h->GetXaxis()->SetTitle( "r.#Delta#phi_{clus-truth}.#sqrt{12}/pitch" );
        h->GetXaxis()->SetRangeUser( -max_det_residual[idet], max_det_residual[idet] );
        h->SetMaximum( h->GetMaximum()*1.2 );

        const int cvid = icut*nLayers[idet] + ilayer + 1;
        cv->cd( cvid );
        h->Draw();

        // fit
        if( h->GetEntries() )
        {
          if( do_fit )
          {
            // const auto result = std::min( Fit( h ), Fit_box( h ) );
            const auto result = Fit_box( h );
            auto f = result._function;
            f->Draw("same");
            auto h = f->GetHistogram();
            auto rms = h->GetRMS();
            auto error = f->GetParError(2);

            Draw::PutText( 0.2, 0.8, Form( "#sigma = %.3g #pm %.3g", rms, error ) );

          } else {

            auto rms = h->GetRMS();
            auto error = h->GetRMSError(2);
            Draw::PutText( 0.2, 0.8, Form( "#sigma = %.3g #pm %.3g", rms, error ) );

          }
        }

        Draw::PutText( 0.2, 0.75, cut_label[icut] );
        Draw::PutText( 0.2, 0.85, Form( "Entries: %.0f", h->GetEntries() ) );

        // draw vertical line at zero
        gPad->Update();
        Draw::VerticalLine( gPad, 0 )->Draw();

      }
    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );

  }

  return pdfFile;

}
