#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/FitUtils.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TH1.h>
#include <TPaveText.h>
#include <TStyle.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

static constexpr double m_jpsi = 3.0969;
static constexpr double m_upsilon = 9.460;
static constexpr double m_e = 0.5110E-3;

static constexpr double m_ref = m_jpsi;

double square( const double x ) { return x*x; }

double fit_function( double* xx, double* par )
{
  const double m = xx[0];
  const double scale = par[0];
  const double r = std::sqrt(1.-square( 2.*m_e/m) );
  const double result = scale*(m/(square(m_ref)-square(m)))*(1+pow(m/m_ref,4))*(std::log((1+r)/(1-r))-r);
  return result;
}

void RadiativeDecay_jpsi(void)
{

  static constexpr double mass_min = 1;
  static constexpr double mass_max = 4;
  static constexpr int nbins = 50*(mass_max - mass_min );
  
  set_style( false );
  const TString tag = "_jpsi_evtgen_acts_full_no_distortion";
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );
  
  const TString pdfFile =  Form( "Figures/RadiativeDecay%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();
    
  {
    // generator level invariant mass
    const TString var( "_track_pairs[]._truth_m" );

    // const TCut pattern_cut;
    const TCut pattern_cut( 
      "_track_pairs[]._tracks[0]._nclusters_mvtx>2"
      "&&_track_pairs[]._tracks[1]._nclusters_mvtx>2"
      "&&_track_pairs[]._tracks[0]._nclusters_tpc>30"
      "&&_track_pairs[]._tracks[1]._nclusters_tpc>30" );
    
    const TString hname( "invMassH_truth" );
    auto h = new TH1F( hname, "", nbins, mass_min, mass_max );
    Utils::TreeToHisto( tree, hname, var, pattern_cut, false );
    h->SetTitle( "" );

    // fit
    auto f = new TF1( "RadDecayF", fit_function, mass_min, mass_max, 1 );
    f->SetNpx( 250*(mass_max - mass_min) );
    f->SetParNames( "scale" );
    f->SetParameter(0, 1 );
    h->Fit( f, "0" );
    f->SetLineColor( 2 );
    
    // draw
    auto cv = new TCanvas( "cv", "cv", 800, 800 );
    cv->SetLeftMargin( 0.15 );
    h->SetMarkerStyle( 20 );
    h->GetXaxis()->SetTitle( "M_{e+e-} (GeV/#it{c}^{2})" );
    h->SetMaximum( 1.2*h->GetMaximum() );
    // h->Draw( "E" );
    h->Draw( "" );
    f->Draw("same");
    gPad->SetLogy( true );
    
    pdfDocument.Add( cv );
    
  }
  
}
