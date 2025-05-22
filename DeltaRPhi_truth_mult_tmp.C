#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <g4eval/TrackingEvaluator_hp.h>

#include <TTree.h>

#include <algorithm>

R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

// tpc layers
static constexpr unsigned int firstlayer_tpc = 7;
static constexpr unsigned int nlayers_tpc = 48;

//_________________________________________________________________________
template< class T>
constexpr T square( const T& x ) { return x*x; }

//_________________________________________________________________________
template< class T>
T delta_phi( const T& phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//_________________________________________________________________________
int count_tpc_clusters( TrackingEvaluator_hp::Container* container )
{
  // loop over tracks
  int out = 0;
  for( const auto& track:container->tracks() )
  {

    // loop over clusters
    out = std::accumulate(
      track._clusters.begin(), track._clusters.end(), out,
      []( int value, const ClusterStruct& cluster )
    { return (cluster._layer >= firstlayer_tpc && cluster._layer < firstlayer_tpc + nlayers_tpc ) ? value+1:value; } );
  }

  return out;
}

//_________________________________________________________________________
TString DeltaRPhi_truth_mult_tmp( TString tag = TString() )
{

  set_style( false );
  gStyle->SetPadLeftMargin(0.13);

  if( tag.IsNull() ) tag = "_hijing_christof" ;

  // select layer (last layer in the TPC)
  static constexpr int layer = 54;
  static constexpr double max_residual = 6;
//   static constexpr int layer = 7;
//   static constexpr double max_residual = 2;

  const TString pdfFile = Form( "Figures/DeltaRPhi_truth_mult%s_%i.pdf", tag.Data(), layer );
  PdfDocument pdfDocument( pdfFile );

  // total tpc clusters per event
  static constexpr float maxmult = 7e5;
  auto h = new TH1F( "h", "h", 100, 0, maxmult );
  h->GetXaxis()->SetTitle( "N_{cluster, TPC}" );
  h->GetXaxis()->SetMaxDigits(3);

  // tpc clusters in tracks per event
  auto htrk = new TH1F( "htrk", "htrk", 100, 0, 7e4 );
  htrk->GetXaxis()->SetTitle( "N_{cluster, TPC, trk}" );
  htrk->GetXaxis()->SetMaxDigits(3);

  // correlate clusters vs clusters in tracks
  auto hcor = new TH2F( "hcor", "hcor", 100, 0, maxmult, 100, 0, 7e4 );

  // residual vs multiplicity
  static constexpr int nmultbins = 10;
  auto profile = new TProfile( "p", "p", nmultbins, 0, maxmult, "s" );
  auto h2 = new TH2F( "h2", "h2", nmultbins, 0, maxmult, 100, -max_residual, max_residual );
  h2->GetXaxis()->SetTitle( "N_{cluster, TPC}" );
  h2->GetXaxis()->SetMaxDigits(3);
  h2->GetYaxis()->SetTitle( "r.#Delta#phi_{track-truth} (cm)" );

  // container
  auto container = new  TrackingEvaluator_hp::Container();

  for( int ifile = 0; ifile<1000; ++ifile )
  {

    // read tpc multiplicity
    int nclusters_tpc;

    {
      TString filename = Form( "DST/CONDOR_Hijing_Christof_simeval/dst_simeval_HijMBPu100_Mar20_1_%04i.root", 1000+ifile );
      std::unique_ptr<TFile> file( TFile::Open( filename ) );
      if( !file ) continue;

      std::unique_ptr<TTree> tree( static_cast<TTree*>( file->Get( "T" ) ) );
      if( !tree ) continue;

      tree->SetBranchAddress( "DST#EVAL#TrackingEvaluator_hp::Container", &container );
      tree->GetEntry(0);

      nclusters_tpc = container->events()[0]._nclusters_tpc;
      h->Fill( nclusters_tpc );
      // std::cout << "DeltaRPhi_truth_mult_tmp - index: " << ifile << " nclusters_tpc: " << nclusters_tpc << std::endl;
    }

    {
      // read tpc track multplicity
      TString filename = Form( "DST/CONDOR_Hijing_Christof/dst_eval_HijMBPu100_Mar20_1_%04i.root", 1000+ifile );

      std::unique_ptr<TFile> file( TFile::Open( filename ) );
      if( !file ) continue;

      std::unique_ptr<TTree> tree( static_cast<TTree*>( file->Get( "T" ) ) );
      if( !tree ) continue;

      tree->SetBranchAddress( "DST#EVAL#TrackingEvaluator_hp::Container", &container );
      tree->GetEntry(0);

      const int nclusters_trk_tpc = count_tpc_clusters( container );
      htrk->Fill( nclusters_trk_tpc );
      hcor->Fill( nclusters_tpc, nclusters_trk_tpc );
      // std::cout << "DeltaRPhi_truth_mult_tmp - index: " << ifile << " nclusters_trk_tpc: " << nclusters_trk_tpc << std::endl;

      // loop over tracks
      for( const auto& track:container->tracks() )
      {

        // embed cut
        if( track._embed != 0 ) continue;

        // momentum cut
        if( track._pt <= 0.5 ) continue;

        // loop over clusters
        for( const auto& cluster:track._clusters )
        {
          // layer cut
          if( cluster._layer != layer ) continue;

          // residual
          const double drphi = cluster._trk_r*delta_phi( cluster._trk_phi - cluster._truth_phi );

          // fill
          h2->Fill( nclusters_tpc, drphi );
          profile->Fill( nclusters_tpc, drphi );
        }
      }

    }

  }

  // plots
  {
    // cluster multiplicity
    auto cv( new TCanvas( "cv", "cv", 800, 800 ) );
    cv->Divide( 1, 2 );
    cv->cd(1);
    gPad->SetRightMargin( 0.12 );
    h->Draw();

    cv->cd(2);
    gPad->SetRightMargin( 0.12 );
    htrk->Draw();
    pdfDocument.Add( cv );
  }

  {
    auto cv( new TCanvas( "cv", "cv", 800, 800 ) );
    cv->SetRightMargin( 0.12 );
    hcor->Draw();
    pdfDocument.Add( cv );
  }

  {
    auto cv( new TCanvas( "cv2", "cv2", 800, 800 ) );
    cv->SetRightMargin( 0.12 );
    h2->Draw();

    profile->SetMarkerColor(2);
    profile->SetLineColor(2);
    profile->Draw( "same" );

    pdfDocument.Add( cv );
  }

  // create TGraph to store resolution vs layer
  auto tg = new TGraphErrors();
  tg->SetName( "residuals" );

  {
    auto cv( new TCanvas( "cv3", "cv3", 800, 800 ) );
    Draw::DivideCanvas( cv, nmultbins, false );
    // fit slices
    for( int i = 0; i < 8; ++i )
      // for( int i = 0; i < nmultbins; ++i )
    {
      const TString hname = Form( "hslice_%i", i );
      auto h = h2->ProjectionY( hname, i+1, i+1 );
      cv->cd(i+1);
      h->Draw();

      const auto result = std::min( Fit( h ), Fit_box( h ) );
      if( result._valid )
      {
        auto f = result._function;
        f->Draw("same");
        auto h = f->GetHistogram();
        auto rms = h->GetRMS();
        auto error = f->GetParError(2);

        Draw::PutText( 0.2, 0.8, Form( "#sigma = %.3g #pm %.3g #mum", rms*1e4, error*1e4 ) );

        tg->SetPoint( i, h2->GetXaxis()->GetBinCenter( i+1 ), rms*1e4 );
        tg->SetPointError( i, 0, error*1e4 );
      }

      // draw vertical line at zero
      gPad->Update();
      Draw::VerticalLine( gPad, 0 )->Draw();

    }
    pdfDocument.Add( cv );
  }

  {
    auto cv( new TCanvas( "cv4", "cv4", 800, 800 ) );
    cv->SetLeftMargin( 0.16 );
    cv->SetRightMargin( 0.12 );

    auto h = new TH1F( "dummy", "", 100, 0, maxmult );
    h->SetMaximum(1.2*Utils::GetMaximum( tg ));
    h->GetXaxis()->SetTitle( "N_{cluster, TPC}" );
    h->GetXaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (track-truth) (#mum)" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();
    tg->SetMarkerStyle(20);
    tg->SetLineColor(1);
    tg->SetMarkerColor(1);
    tg->Draw("P");

    pdfDocument.Add( cv );
  }

  return pdfFile;

}
