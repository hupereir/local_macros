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

static constexpr float z_ref = 5.;

static constexpr int nxbins = 160;

namespace { 
  
  template<class T> 
    constexpr T square( const T& x ) { return x*x; } 

  template< class T>
    T normalize_phi( T phi )
  {
    while( phi >= 2*M_PI ) phi -= 2*M_PI;
    while( phi < 0 ) phi += 2*M_PI;
    return phi;
  }

  // convert phi, r, z histogram into x,y,z
  TH2* to_cartesian( TH3* h_radial, const TString name )
  {
    
    // get maximum radius
    const auto r_min = 30;
    // const auto r_min = h_radial->GetYaxis()->GetBinLowEdge( 2 );
    const auto r_max = h_radial->GetYaxis()->GetBinUpEdge( h_radial->GetYaxis()->GetNbins()-1 );
    
    // get cartesian
    auto h_cartesian = new TH2F( name, name, nxbins, -r_max, r_max, nxbins, -r_max, r_max );
    
    // loop over bins
    for( int ix = 0; ix < nxbins; ++ix )
      for( int iy = 0; iy < nxbins; ++iy )
    {
      const auto x = h_cartesian->GetXaxis()->GetBinCenter( ix+1 );
      const auto y = h_cartesian->GetYaxis()->GetBinCenter( iy+1 );
      const auto phi = normalize_phi( std::atan2(y, x) );
      const auto r = std::sqrt( square(x) + square(y) );
      if( r < r_min ) continue;
      if( r > r_max ) continue;
      
      
      const auto distortion = h_radial->Interpolate( phi, r, z_ref );
      h_cartesian->SetBinContent( ix+1, iy+1, distortion );
      
    }  
    return h_cartesian;
  }

  double get_sector_phi( int isec ) { return isec*M_PI/6 + M_PI/12; }

  void draw_sector_boundaries()
  {
    
    // draw sector boundaries
    for( int i = 0; i < 12; ++i )
    {
      const auto phi( get_sector_phi(i) + M_PI/12 );
      const auto cphi = std::cos(phi);
      const auto sphi = std::sin(phi);
      static const double r_min = 20;
      static const double r_max = 80;
      
      auto line = new TLine( r_min*cphi, r_min*sphi, r_max*cphi, r_max*sphi );
      line->SetLineStyle( 2 );
      line->SetLineWidth( 1 );
      line->SetLineColor( 1 );
      line->Draw();
    }
    
    // also draw micromegas sector
    {
      constexpr int i = 3;
      const auto phi( get_sector_phi(i) );
      const auto cphi = std::cos(phi);
      const auto sphi = std::sin(phi);
      static const double r_min = 20;
      static const double r_max = 80;
      
      auto line = new TLine( r_min*cphi, r_min*sphi, r_max*cphi, r_max*sphi );
      line->SetLineStyle( 2 );
      line->SetLineWidth( 1 );
      line->SetLineColor( 2 );
      line->Draw();
    }    
  
  }
  
}

//_____________________________________________________________________________
TString DrawDistortionMap_xy()
{

  set_style( false );
  
  gStyle->SetPalette( kBird );
  // gStyle->SetPalette( kBlueGreenYellow );

  // open TFile
//   const TString tag = "_fluct_average-coarse";
//   const TString inputfile = "distortion_maps/fluct_average-coarse.root";

//   const TString tag = "_full_realistic_micromegas_all-coarse";
//   const TString inputfile = "distortion_maps_rec/Distortions_full_realistic_micromegas_all-coarse.root";
  
//   const TString tag = "_full_realistic_micromegas_all-coarse";
//   const TString inputfile = "distortion_maps_rec/Distortions_full_realistic_micromegas_all-coarse.root";
  
  const TString tag = "_full_realistic_micromegas_all-coarse";
  const TString inputfile = "distortion_maps_rec/Distortions_full_realistic_micromegas_all-coarse.root";

//   const TString tag = "full_realistic_micromegas_mm-coarse_extrapolated";
//   const TString inputfile = "distortion_maps_rec/Distortions_full_realistic_micromegas_mm-coarse_extrapolated.root";

//   const TString tag = "full_realistic_micromegas_truth-coarse";
//   const TString inputfile = "distortion_maps_rec/Distortions_full_realistic_micromegas_truth-coarse.root";

  auto f = TFile::Open( inputfile );
  if( !f ) return TString();

  TString pdfFile( Form( "Figures/DistortionMap_xy%s.pdf", tag.Data() ) );
  PdfDocument pdfDocument( pdfFile );
  
  using histogram_struct_t = std::tuple<TH3*, TString, TString>;
//   std::array<histogram_struct_t, 3> histogram_structs = {{
//     { dynamic_cast<TH3*>(f->Get("hDistortionP_rec")), "hIntDistortionR_xy", "r#Delta#phi (cm)" },
//     { dynamic_cast<TH3*>(f->Get("hDistortionR_rec")), "hIntDistortionP_xy", "#Deltar (cm)" },
//     { dynamic_cast<TH3*>(f->Get("hDistortionZ_rec")), "hIntDistortionZ_xy", "#Deltaz (cm)" }
//   }};

  std::array<histogram_struct_t, 3> histogram_structs = {{
    { dynamic_cast<TH3*>(f->Get("hIntDistortionP")), "hIntDistortionR_xy", "r#Delta#phi (cm)" },
    { dynamic_cast<TH3*>(f->Get("hIntDistortionR")), "hIntDistortionP_xy", "#Deltar (cm)" },
    { dynamic_cast<TH3*>(f->Get("hIntDistortionZ")), "hIntDistortionZ_xy", "#Deltaz (cm)" }
  }};
  
  TCanvas* cv( new TCanvas( "cv", "cv", 1000, 300 ) );
  cv->Divide( 3, 1 );
  int cvid = 0;
  for( const auto& [h3, name, label]:histogram_structs )
  {
    auto h = to_cartesian( h3, name );
    h->GetXaxis()->SetTitle( "x (cm)" );
    h->GetYaxis()->SetTitle( "y (cm)" );
    h->GetZaxis()->SetTitle( label );
    h->GetZaxis()->SetTitleOffset( 1.5 );
    h->SetTitle("");
    cv->cd( ++cvid );
    h->Draw("colz");
    gPad->SetRightMargin( 0.18 );
    gPad->Update();
    
    draw_sector_boundaries();
  
  }

  pdfDocument.Add( cv );
  return pdfFile;

}
