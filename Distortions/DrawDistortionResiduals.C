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

R__LOAD_LIBRARY(libRootUtilBase.so)

constexpr int m_phibins = 36;
constexpr int m_rbins = 16;
constexpr int m_zbins = 80;

TString DrawDistortionResiduals()
{
  set_style( false );
  
  /* bin definitions */
  /*
  For Me
  x axis: phi, 36 bins 0 to 2pi
  y axis: r, 16 bins, 20 to 78
  z axis: z, 80 bins -105.5, +105.5
  */

  // open TFile
  const TString tag = "_Hijing_Micromegas_50kHz_truth_notpc_all";
  // const TString tag = "_Hijing_Micromegas_50kHz_truth_notpc";
  auto f = TFile::Open( Form( "Rootfiles/Distortions_drphi_full%s.root", tag.Data() ) );
  if( !f ) return TString();

  TString pdfFile( Form( "Figures/DistortionResidual%s.pdf", tag.Data() ) );
  PdfDocument pdfDocument( pdfFile );

  // define bins
  static constexpr int phibin = 12;
  static constexpr int zbin = 40;
  
  {
    auto cv = new TCanvas( "cv", "cv", 800, 800 );
    Draw::DivideCanvas( cv, m_rbins );
    
    for( int i = 0; i < m_rbins; ++i )
    {
      auto hname = Form( "residual_p%i_r%i_z%i", phibin, i, zbin );
      auto h = dynamic_cast<TH1*>( f->Get( hname ) );
      cv->cd( i+1 );
      if( h ) h->Draw();
    }
    
    pdfDocument.Add( cv );
    
  }
  
  return pdfFile;
    
}

