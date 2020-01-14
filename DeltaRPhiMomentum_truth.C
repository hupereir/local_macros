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
TString DeltaRPhiMomentum_truth( TString tag = TString() )
{

  set_style( false );

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  const int layer = 7;

  // initial guess for max residual
  Float_t max_residual = 1.5;

  if( tag.IsNull() ) tag = "_low_momentum_truth_notpc_noouter" ;
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*.root", tag.Data(), tag.Data() );
  const TString pdfFile = Form( "Figures/DeltaRPhiMomentum_truth%s_%i.pdf", tag.Data(), layer );

  std::cout << "DeltaRPhi_truth - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhi_truth - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "_clusters._r*DeltaPhi(_clusters._trk_phi - _clusters._truth_phi)" );
  const TString var2d = Form( "%s:_clusters._truth_pt", var.Data() );
  const TCut layer_cut( Form( "_clusters._layer == %i", layer ) );

  // create TGraph to store resolution vs layer
  auto tg = new TGraphErrors();
  tg->SetName( "residuals" );

  // 2d histogram
  TString hname( "deltarphi" );

  // momentum range
  const double ptmin = 0.5;
  const double ptmax = 5;
  constexpr int nptbins = 40;

  // optimize max esidual
  for( int i=0; i<3; ++i )
  {
    const TString hname( "h" );
    std::unique_ptr<TH1> h1( new TH1F( hname, "", 500, -max_residual, max_residual ) );
    Utils::TreeToHisto( tree, hname, var, layer_cut, false );
    max_residual = 5*h1->GetRMS();
  }

  auto h2d( new TH2F( hname, "", nptbins, ptmin, ptmax, 100, -max_residual, max_residual ) );
  Utils::TreeToHisto( tree, hname, var2d, layer_cut, false );
  std::cout << "Entries: " << h2d->GetEntries() << endl;

  // canvas
  const TString cvName( "cv0" );
  auto cv = new TCanvas( cvName, cvName, 800, 800 );
  Draw::DivideCanvas( cv, nptbins, false );

  // loop over bins and fit
  for( int ibin = 0; ibin < nptbins; ++ibin )
  {
    const auto hname = Form( "h_%i", ibin );
    auto h = h2d->ProjectionY( hname, ibin+1, ibin+1 );
    cv->cd( ibin+1 );
    h->SetMaximum( 1.1*h->GetMaximum() );
    h->Draw();

    // auto result = Fit( h );
    const auto result = std::min( Fit( h ), Fit_box( h ) );
    std::cout << "valid: " << result.valid() << std::endl;
    auto f = result._function;
    f->Draw( "same" );

    {
      auto h = f->GetHistogram();
      auto sigma = h->GetRMS()*1e4;
      auto sigmaError = f->GetParError(2)*1e4;

      Draw::PutText( 0.2, 0.8, Form( "#sigma = %.3g #pm %.3g #mum", sigma, sigmaError ) );
      tg->SetPoint( ibin, h2d->GetXaxis()->GetBinCenter( ibin+1 ), sigma );
      tg->SetPointError( ibin, 0, sigmaError );
    }

  }

  pdfDocument.Add( cv );

  // TGraph
  cv = new TCanvas( "cvtg", "cvtg", 800, 800 );
  cv->SetLeftMargin( 0.17 );

  auto h = new TH1F( "dummy", "", 100, ptmin, ptmax );
  h->SetMinimum(0);
  h->SetMaximum(max_residual*1e4/2);
  h->GetXaxis()->SetTitle( "#it{p}_{T} (GeV/#it{c})" );
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
