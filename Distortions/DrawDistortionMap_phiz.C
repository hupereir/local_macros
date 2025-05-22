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

//_________________________________________________
TH2* mask_empty_bins( TH2* ref, TH2* mask )
{
  const auto minimum = ref->GetMinimum();
  const auto maximum = ref->GetMaximum();
  if( minimum > 0 || maximum < 0 ) return ref;

  auto copy = static_cast<TH2*>( ref->Clone() );
  const auto minimum_scaled=minimum*1.2;
  for( int ix = 0; ix < ref->GetNbinsX(); ++ix )
    for( int iy = 0; iy < ref->GetNbinsY(); ++iy )
  {
    auto masked = mask->GetBinContent( ix+1, iy+1 );
    if( masked > 0 ) { copy->SetBinContent( ix+1, iy+1, ref->GetBinContent(ix+1, iy+1) ); }
    else copy->SetBinContent( ix+1, iy+1, minimum_scaled );
  }
  copy->GetZaxis()->SetRangeUser(minimum, maximum);
  return copy;
}

//_________________________________________________
TString DrawDistortionMap_phiz()
{

  set_style( false );
  gStyle->SetOptStat(0);

//   const TString tag = "_average_minus_static_projected";
//   const TString inputfile = "distortion_maps/average_minus_static_distortion_projected.root";
//   const TString pdfFile( Form( "Figures/DistortionMap%s_phiz.pdf", tag.Data() ) );
//   const bool use_phi_as_radian = false;

  const TString tag = "_average_minus_static_extrapolated";
  const TString inputfile = "distortion_maps/average_minus_static_distortion_extrapolated.root";
  const TString pdfFile( Form( "Figures/DistortionMap%s_phiz.pdf", tag.Data() ) );
  const bool use_phi_as_radian = false;

//   const TString tag = "_average_minus_static_converted";
//   const TString inputfile = "distortion_maps/average_minus_static_distortion_converted.root";
//   const TString pdfFile( Form( "Figures/DistortionMap%s_phiz.pdf", tag.Data() ) );
//   const bool use_phi_as_radian = false;

  std::cout << "DrawDistortionMap - inputfile: " << inputfile << std::endl;
  std::cout << "DrawDistortionMap - pdfFile: " << pdfFile << std::endl;

  // open TFile
  auto f = TFile::Open( inputfile );
  if( !f ) return TString();

  PdfDocument pdfDocument( pdfFile );

  const TString prefix[] = {
    "hIntDistortionR",
    "hIntDistortionP",
    "hIntDistortionZ"
  };

  const TString suffix[] = {
    "_negz",
    "_posz"
  };

  const TString label[] = {
    "#Deltar (cm)",
    "r#Delta#phi (cm)",
    "#Deltaz (cm)"
  };

  const TString cvname = "cv";
  TCanvas* cv( new TCanvas( cvname, cvname, 1200, 800 ) );
  cv->Divide( 3, 2 );

  auto select_range = []( TAxis* a, double x_min, double x_max )
  {
    int bin_min = a->FindBin( x_min );
    int bin_max = a->FindBin( x_max );
    a->SetRangeUser(x_min, x_max);
    return 1./(bin_max - bin_min+1);
  };

  for( int isuffix = 0; isuffix < 2; ++isuffix )
  {

    for( int ih = 0; ih < 3; ++ih )
    {
      const auto hname = Form( "%s%s", prefix[ih].Data(), suffix[isuffix].Data() );
      const auto maskname = Form( "hentries%s", suffix[isuffix].Data() );

      auto h3= dynamic_cast<TH3*>(f->Get(hname));
      auto hmask3 = dynamic_cast<TH3*>(f->Get(maskname));

      const double r_ref = h3->GetYaxis()->GetBinCenter( h3->GetYaxis()->FindBin( 40 ) );
      const double scale = select_range( h3->GetYaxis(), r_ref, r_ref );

      TH2* h = static_cast<TH2*>(h3->Project3D( "zx" ));
      h->SetTitle( "" );
      h->Scale( scale );
      h->GetXaxis()->SetTitle( "#phi (rad)" );
      h->GetYaxis()->SetTitle( "z (cm)" );
      h->GetZaxis()->SetTitle( label[ih] );
      h->GetZaxis()->SetTitleOffset( 1.8 );

      if( hmask3 )
      {
        select_range( hmask3->GetYaxis(), r_ref, r_ref );
        auto hmask2 = static_cast<TH2*>( hmask3->Project3D( "zx" ) );
        h = mask_empty_bins( h, hmask2 );
      }

      // rescale as rdphi
      if( use_phi_as_radian && ih ==1 )
      {
        for( int iy = 0; iy < h->GetNbinsY(); ++iy )
        {
          for( int ix = 0; ix < h->GetNbinsX(); ++ix )
          { h->SetBinContent( ix+1, iy+1, h->GetBinContent( ix+1, iy+1 )*r_ref ); }
        }
      }

      cv->cd(ih+3*isuffix+1);
      gPad->SetRightMargin( 0.22 );

      h->Draw("colz");
    }
  }

  pdfDocument.Add( cv );
  return pdfFile;

}
