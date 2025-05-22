#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TCut.h>
#include <TH1.h>
#include <TH2.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

bool exclude_points = false;

//______ _____________________________________________________________________
double fit_function( double* x, double* par )
{
  const double z = x[0];
  const double alpha = par[0];
  if( exclude_points && std::abs(z)<5 )
  {
    TF1::RejectPoint(true);
    return 0;
  }

  if( z > 0 ) return alpha*(105.5-z);
  else return alpha*(-105.5-z);
}

//____________________________________________________________________________
TString DeltaZ_z()
{

  set_style( false );

  // input files
//   const TString tag = "_full_realistic-scaled2-new";
//   const TString subtag = "";
//   const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );
//   const TString tag = "_full_realistic-scaled2-new";


  // input files
  // const TString tag = "_realistic_acts_truth_notpc_nodistortion";
  const TString tag = "_realistic_acts_truth_notpc_vdrift2";
  const TString inputFile = Form( "DST/CONDOR%s/dst_reco*_??.root", tag.Data() );
  const TString pdfFile = Form( "Figures/DeltaZ_z%s.pdf", tag.Data());

  std::cout << "DeltaZ - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaZ - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  bool do_fit = false;

  // variable names
  const TString var( "_tracks._clusters._trk_z - _tracks._clusters._z" );
  const TString var2d = Form( "%s:_tracks._clusters._z", var.Data() );
  const TCut momentum_cut = ("_tracks._pt>0.5" );
  const TCut pattern_cut(
    "_tracks._truth_pt>0.5"
    "&&_tracks._nclusters_mvtx>0"
    "&&_tracks._nclusters_intt>=2"
    "&&_tracks._nclusters_micromegas==2"
    );

//   const TCut radius_cut = ( "_tracks._clusters._r < 50" );
  const TCut radius_cut;
  const TCut dz_cut( "fabs(_tracks._clusters._trk_z - _tracks._clusters._z)<10" );

  // select only clusters in the TPC
  const TCut detector_cut( "_tracks._clusters[]._layer>=7 &&_tracks._clusters[]._layer<55 " );

  const TString hname( "hdeltaz_z" );
  static constexpr float max_dz = 0.3;
  // static constexpr float max_dz = 3;
  auto h2d( new TH2F( hname, "", 100, -105.5, 105.5, 100, -max_dz, max_dz ) );
  h2d->GetXaxis()->SetTitle( "z (cm)" );
  h2d->GetYaxis()->SetTitle( "#Deltaz (track - cluster) (cm)" );
  Utils::TreeToHisto( tree, hname, var2d, momentum_cut&&pattern_cut&&detector_cut&&dz_cut&&radius_cut, false );

  const TString pname( "pdeltaz_z" );
  auto p( new TProfile( pname, "", 100, -105.5, 105.5 ) );
  Utils::TreeToHisto( tree, pname, var2d, momentum_cut&&pattern_cut&&detector_cut&&dz_cut&&radius_cut, false );

  // canvas
  std::unique_ptr<TCanvas> cv( new TCanvas( "cvtg", "cvtg", 800, 600 ) );
  cv->SetRightMargin(.2);
  h2d->Draw( "colz");

  // auto p = h2d->ProfileX();
  p->SetLineColor(2);
  p->SetMarkerColor(2);
  p->Draw( "same" );

  Draw::HorizontalLine( cv.get(), 0 )->Draw();

  if( do_fit )
  {
    const TString fname = "fit_function";
    auto f = new TF1( fname, ::fit_function, -80, 80, 1 );
    // auto f = new TF1( fname, ::fit_function, -105.5, 105.5, 1 );
    f->SetParameter( 0, 1.0 );

    exclude_points = true;
    // p->Fit( f, "0" );
    h2d->Fit( f, "0" );

    exclude_points = false;
    f->SetLineColor( 2 );
    f->Draw("same");

    auto legend = new TLegend( 0.5, 0.85, 0.78, 0.9, "", "NDC" );
    //   legend->SetFillColor(0);
    //   legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry( f, Form( "#alpha = %.3g #pm %.2g", f->GetParameter(0), f->GetParError(0) ), "" );
    legend->Draw();
  }

  pdfDocument.Add( cv.get() );

  return pdfFile;

}
