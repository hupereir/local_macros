#include <RootUtil/Draw.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH3.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

// binning and range
constexpr int m_phibins = 36;
constexpr int m_rbins = 16;
constexpr int m_zbins = 80;

// phi range
constexpr float m_phimin = 0;
constexpr float m_phimax = 2.*M_PI;

// r range
constexpr float m_rmin = 20;
constexpr float m_rmax = 78;

// z range
constexpr float m_zmin = -105.5;
constexpr float m_zmax = 105.5;

static constexpr int isec_rec = 3;
static constexpr double phi_rec = M_PI*isec_rec/6 + M_PI/12;


//__________________________________________________________________________
double get_phi( double phi )
{ return phi >= 0 ? phi: (phi+2.*M_PI); }

//__________________________________________________________________________
TString ScanClusters_coarse()
{
  set_style( false );
  // const TString tag = "_flat_truth_micromegas_nominal";
  // const TString tag = "_flat_truth_micromegas_corrected_mm-coarse_extrapolated-new";
  const TString tag = "_flat_truth_micromegas_corrected_mm-coarse_extrapolated-new2";

  // const TString inputFile = Form( "DST/CONDOR%s/dst_eval*_1?.root", tag.Data() );
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval*.root", tag.Data() );
  const TString pdfFile = Form( "Figures/ScanClusters_coarse%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // variable
  const TString var( "_clusters._z:_clusters._r:get_phi(_clusters._phi)" );

  // load input
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // cluster scan
  auto h = new TH3F( "scan", "scan",
    m_phibins, m_phimin, m_phimax,
    m_rbins, m_rmin, m_rmax,
    m_zbins, m_zmin, m_zmax );

  Utils::TreeToHisto( tree, h->GetName(), var, TCut(), false );
  h->SetTitle( "" );
  h->GetXaxis()->SetTitle( "#phi (rad)" );
  h->GetYaxis()->SetTitle( "r (cm)" );
  h->GetZaxis()->SetTitle( "z (cm)" );

  // project on layers
  for( int ir = 0; ir < m_rbins; ++ir )
  {
    const auto hname = Form( "%s_%i", h->GetName(), ir );

    // layer range
    h->GetYaxis()->SetRange( ir+1, ir+1 );
    const auto hlayer( h->Project3D( "zx" ) );
    hlayer->SetName( hname );
    hlayer->SetTitle("");

    const auto cvname = Form( "cv_%i", ir );
    auto cv = new TCanvas( cvname, cvname, 800, 800 );
    cv->SetLeftMargin( 0.14 );
    cv->SetRightMargin( 0.20 );

    hlayer->Draw( "colz" );
    hlayer->GetYaxis()->SetTitleOffset( 1.35 );

    Draw::PutText( 0.2, 0.8, Form( "r = %.2fcm - entries = %.0f",
      h->GetYaxis()->GetBinCenter( ir+1 ),
      hlayer->GetEntries() ) );

    // draw sector boundaries
    for( int i = 0; i < 12; ++i )
    {
      const double phi = (i+1)*M_PI/6;
      Draw::VerticalLine( cv, phi )->Draw();
    }

    {
      auto line = Draw::VerticalLine( cv, phi_rec );
      line->SetLineColor( 2 );
      line->Draw();
    }

    pdfDocument.Add( cv );

  }

  return pdfFile;

}
