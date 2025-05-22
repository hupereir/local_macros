#include <RootUtil/Draw.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>


#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <cmath>

R__LOAD_LIBRARY(libRootUtilBase.so)

//_______________________________________________________________
float square( float value ) { return value*value; }

//_______________________________________________________________
float get_r( float x, float y ) { return std::sqrt( square(x)+square(y) ); }

//____________________________________________________________________________
TString make_run_label( const std::vector<int>& runlist )
{
  if( runlist.empty() ) return TString();
  if( runlist.size() == 1 ) return Form( "run %i", runlist[0] );
  return Form( "runs %i-%i",
    *std::min_element( runlist.begin(), runlist.end()),
    *std::max_element( runlist.begin(), runlist.end()) );
}

//____________________________________________________________________________
TString make_run_postfix( const std::vector<int>& runlist )
{
  if( runlist.empty() ) return TString();
  if( runlist.size() == 1 ) return Form( "_%08i", runlist[0] );
  return Form( "_%08i-%08i",
    *std::min_element( runlist.begin(), runlist.end()),
    *std::max_element( runlist.begin(), runlist.end()) );
}

//_______________________________________________________________
double fit_function( double* x, double* par )
{
  const double z = x[0];
  const double dv = par[0];
  const double dz_minus = par[1];
  const double dz_plus = par[2];

  static const double zmin = -110;
  static const double zmax = +110;


  return (z)<0 ?
    dz_minus + dv*z:
    dz_plus + dv*z;
}

//_______________________________________________________________
void TpotMatchingDeltaZ_fit()
{

  set_style( false );

  const TString tag = "_CombinedDataReconstruction_corrected-2";
  TString postfix;
  TString run_label;

  // const std::vector<int> runlist = { 53199 };
  const std::vector<int> runlist = { 53285 };
  run_label = make_run_label( runlist );
  postfix = make_run_postfix( runlist );

  const TString rootFilename = Form( "Rootfiles/TpotMatchingDeltaZ%s%s.root", tag.Data(),postfix.Data());
  auto tfile = TFile::Open(rootFilename);

  const TString pdfFile = Form( "Figures/TpotMatchingDeltaZ_fit%s%s.pdf", tag.Data(),postfix.Data());
  PdfDocument pdfDocument( pdfFile );

  auto h = static_cast<TH2*>( tfile->Get("h_dz") );

  auto p = h->ProfileX("p");
  p->SetMarkerStyle(20);

  auto h_mean =  static_cast<TH1*>( tfile->Get("h_dz_1") );

  auto f = new TF1( "fit_function", fit_function, -50, 50, 3 );
  // auto f = new TF1( "fit_function", fit_function, 10, 100, 2 );
  f->SetLineColor(2);
  f->SetLineWidth(2);
  // f->SetParameters(0,0,0);
  f->SetParameters(0.04,-3, 3 );

  h->Fit( f, "0R" );
  // h_mean->Fit( f, "0R" );

  {

    auto cv = new TCanvas( "cv", "cv", 800, 800 );
    cv->cd(1);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);

    h->Draw("colz" );
    gPad->Update();
    Draw::PutText(0.15, 0.92, Form( "%s, entries: %.0f", run_label.Data(), h->GetEntries()));
    Draw::HorizontalLine(gPad, 0)->Draw();

    p->Draw("same");
    h_mean->Draw("same");
    f->Draw("same");

    pdfDocument.Add(cv);
  }
}
