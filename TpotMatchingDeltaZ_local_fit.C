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
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>
#include <micromegas/MicromegasDefs.h>

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
  const double dz = par[1];
  return dz + dv*z;
}

//_______________________________________________________________
void TpotMatchingDeltaZ_local_fit()
{
  MicromegasMapping mapping;

  set_style( false );

  const TString tag = "_CombinedDataReconstruction_corrected-2";
  TString postfix;
  TString run_label;

  // const std::vector<int> runlist = { 53199 };
  const std::vector<int> runlist = { 53285 };
  run_label = make_run_label( runlist );
  postfix = make_run_postfix( runlist );

  const TString rootFilename = Form( "Rootfiles/TpotMatchingDeltaZ_local%s%s.root", tag.Data(),postfix.Data());
  auto tfile = TFile::Open(rootFilename);

  const TString pdfFile = Form( "Figures/TpotMatchingDeltaZ_local_fit%s%s.pdf", tag.Data(),postfix.Data());
  PdfDocument pdfDocument( pdfFile );

  // local coordinates
  auto cv = new TCanvas( "cv2", "cv2", 800, 600 );
  cv->Divide(4, 2);

  // loop over tiles
  for( int itile = 0; itile < MicromegasDefs::m_ntiles; ++itile )
  {

    const auto hname = Form( "h_z_%i", itile );
    auto h = static_cast<TH2*>( tfile->Get(hname) );

    const auto fname = Form( "fit_function_%i", itile );
    auto f = new TF1( fname, fit_function, -20, 20, 2 );
    f->SetLineColor(2);
    f->SetLineWidth(2);
    f->SetParameters(0.04,0);

    h->Fit( f, "0R" );

    cv->cd(itile+1);
    h->Draw("col");
    f->Draw("same");

    // get detector name
    const auto hitsetkey = MicromegasDefs::genHitSetKey(56, MicromegasDefs::SegmentationType::SEGMENTATION_Z, itile);
    const auto detname = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

    gPad->SetTopMargin(0.1);
    Draw::PutText( 0.15, 0.92, Form( "%s - %s, entries: %.0f", detname.c_str(), run_label.Data(), h->GetEntries() ));

    pdfDocument.Add(cv);
  }
}
