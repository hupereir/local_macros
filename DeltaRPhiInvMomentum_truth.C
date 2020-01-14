#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TDirectory.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
Float_t DeltaPhi( Float_t phi )
{
  if( phi >= 2*M_PI ) return phi - 2*M_PI;
  else if( phi <= -2*M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
TString DeltaRPhiInvMomentum_truth( TString tag = TString() )
{

  set_style( kFALSE );

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  const std::array<Float_t, nDetectors> maxDetResidual = { 0.003, 0.01, 1.2, 1.2, 1.2, 1.5 };
  const int layer = 7;
  const auto maxResidual = 0.5;

  // pdf output
  if( tag.IsNull() ) tag = "_5k_truth_low_momentum_notpc_noouter" ;
  const TString inputFile = Form( "DST_afs/dst_eval%s.root", tag.Data() );
  const TString pdfFile = Form( "Figures/DeltaRPhiInvMomentum_truth%s_%i.pdf", tag.Data(), layer );

  std::cout << "DeltaRPhi_truth - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhi_truth - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // variable names
  const TString var( "_clusters._trk_r*DeltaPhi(_clusters._trk_phi - _clusters._truth_phi)" );
  const TString var2d = Form( "%s:1.0/_clusters._truth_pt", var.Data() );
  const TCut layerCut( Form( "_clusters._layer == %i", layer ) );

  // create TGraph to store resolution vs layer
  auto tg = new TGraphErrors();
  tg->SetName( "residuals" );

  // 2d histogram
  TString hname( "deltarphi" );
  const double ptmin = 0.5;
  const double ptmax = 3;
  const int nptbins = 12;
  auto h2d( new TH2F( hname, "", nptbins, 1./ptmax, 1./ptmin, 50, -maxResidual, maxResidual ) );
  Utils::TreeToHisto( tree, hname, var2d, layerCut, kFALSE );
  std::cout << "Entries: " << h2d->GetEntries() << endl;

  // canvas
  const TString cvName( "cv0" );
  auto cv = new TCanvas( cvName, cvName, 800, 800 );
  Draw::DivideCanvas( cv, nptbins, kFALSE );

  // loop over bins and fit
  for( int ibin = 0; ibin < nptbins; ++ibin )
  {
    const auto hname = Form( "h_%i", ibin );
    TH1* h = h2d->ProjectionY( hname, ibin+1, ibin+1 );
    cv->cd( ibin+1 );
    h->Draw();

    auto f = Fit( h );
    Draw::PutText( 0.2, 0.8, Form( "#sigma = %.3g #pm %.3g #mum", f->GetParameter(2)*1e4, f->GetParError(2)*1e4 ) );
    tg->SetPoint( ibin, h2d->GetXaxis()->GetBinCenter( ibin+1 ), f->GetParameter(2)*1e4 );
    tg->SetPointError( ibin, 0, f->GetParError(2)*1e4 );

  }

  pdfDocument.Add( cv );

  // TGraph
  cv = new TCanvas( "cvtg", "cvtg", 800, 800 );
  cv->SetLeftMargin( 0.17 );

  auto h = new TH1F( "dummy", "", 100, 1./ptmax, 1./ptmin );
  h->SetMinimum(0);
  h->SetMaximum(maxResidual*1e4 / 2);
  h->GetXaxis()->SetTitle( "1/#it{p}_{T} ((GeV/#it{c})^{-1})" );
  h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (track-truth) (#mum)" );
  h->GetYaxis()->SetTitleOffset( 1.7 );
  h->Draw();

  tg->SetMarkerStyle(20);
  tg->SetLineColor(1);
  tg->SetMarkerColor(1);
  tg->Draw("P");

  pdfDocument.Add( cv );
  return pdfFile;

}
