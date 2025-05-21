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
TString DrawDistortionMap()
{

  set_style( false );
  gStyle->SetOptStat(0);

//   const TString tag = "_static_only";
//   const TString inputfile = "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.7/share/calibrations/distortion_maps/static_only_inverted_10-new.root";
//   const TString pdfFile( Form( "Figures/DistortionMap%s.pdf", tag.Data() ) );
//   const bool use_phi_as_radian = false;

  const TString tag = "_average";
  const TString inputfile = "/sphenix/tg/tg01/jets/bkimelman/BenProduction/Feb21_2025/Laminations_run2pp_ana464_2024p011_v001-00053534.root";
  const TString pdfFile( Form( "Figures/DistortionMap%s.pdf", tag.Data() ) );
  const bool use_phi_as_radian = false;

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

  auto select_range = []( TAxis* a, double x_min, double x_max )
  {
    int bin_min = a->FindBin( x_min );
    int bin_max = a->FindBin( x_max );
    a->SetRangeUser(x_min, x_max);
    return 1./(bin_max - bin_min+1);
  };

  if( true )
  {

    const TString cvname = "cv";
    TCanvas* cv( new TCanvas( cvname, cvname, 1200, 800 ) );
    cv->Divide( 3, 2 );

    // 2d plots
    for( int isuffix = 0; isuffix < 2; ++isuffix )
    {

      for( int ih = 0; ih < 3; ++ih )
      {
        const auto hname = Form( "%s%s", prefix[ih].Data(), suffix[isuffix].Data() );
        const auto maskname = Form( "hentries%s", suffix[isuffix].Data() );

        // deltaR vs R and z
        auto h3= dynamic_cast<TH3*>(f->Get(hname));
        auto hmask3 = dynamic_cast<TH3*>(f->Get(maskname));

        h3->GetYaxis()->SetRangeUser( 30, h3->GetYaxis()->GetXmax() );
        // const double scale = select_range( h3->GetXaxis(), -1.742 + 2*M_PI,-1.43979+2*M_PI);
        const double scale = select_range( h3->GetXaxis(), -2.27002+2*M_PI,-1.9673+2*M_PI);

        auto h = static_cast<TH2*>(h3->Project3D( "yz" ));
        h->SetTitle( "" );

        h->Scale( scale );
        h->GetXaxis()->SetTitle( "z (cm)" );
        h->GetYaxis()->SetTitle( "r (cm)" );
        h->GetZaxis()->SetTitle( label[ih] );
        h->GetZaxis()->SetTitleOffset( 1.8 );

        if( hmask3 )
        {
          hmask3->GetYaxis()->SetRangeUser( 30, hmask3->GetYaxis()->GetXmax() );
          // const double scale = select_range( hmask3->GetXaxis(), -1.742 + 2*M_PI,-1.43979+2*M_PI);
          const double scale = select_range( hmask3->GetXaxis(), -2.27002+2*M_PI,-1.9673+2*M_PI);
          auto hmask2 = static_cast<TH2*>( hmask3->Project3D( "yz" ) );
          h = mask_empty_bins( h, hmask2 );
        }

        // rescale as rdphi
        if( use_phi_as_radian && ih ==1 )
        {
          for( int iy = 0; iy < h->GetNbinsY(); ++iy )
          {
            const double r = h->GetYaxis()->GetBinCenter( iy+1 );
            for( int ix = 0; ix < h->GetNbinsX(); ++ix )
            { h->SetBinContent( ix+1, iy+1, h->GetBinContent( ix+1, iy+1 )*r ); }
          }
        }

        cv->cd(ih+3*isuffix+1);
        gPad->SetRightMargin( 0.22 );

        h->Draw("colz");
        Draw::HorizontalLine( gPad, 30 )->Draw();
      }
    }

    pdfDocument.Add( cv );
  }

  if( true )
  {
    const TString cvname = "cv1";
    TCanvas* cv( new TCanvas( cvname, cvname, 1200, 800 ) );
    cv->Divide( 3, 2 );

    // 1d plots
    for( int isuffix = 0; isuffix < 2; ++isuffix )
    {

      for( int ih = 0; ih < 3; ++ih )
      {
        const auto hname = Form( "%s%s", prefix[ih].Data(), suffix[isuffix].Data() );
        const auto maskname = Form( "hentries%s", suffix[isuffix].Data() );

        // deltaR vs R and z
        auto h3= dynamic_cast<TH3*>(f->Get(hname));
        h3->GetYaxis()->SetRangeUser( 30, h3->GetYaxis()->GetXmax() );
        const double scale_phi = select_range( h3->GetXaxis(), -2.27002+2*M_PI,-1.9673+2*M_PI);
//         const double scale_z = isuffix ?
//           select_range( h3->GetZaxis(), 0, 100):
//           select_range( h3->GetZaxis(), -100, 0);

        const double scale_z = isuffix ?
          select_range( h3->GetZaxis(), 50, 52):
          select_range( h3->GetZaxis(), -52, -50);

        auto h = static_cast<TH1*>(h3->Project3D( "y" ));
        h->SetTitle( "" );

        h->Scale( scale_phi*scale_z );
        h->GetXaxis()->SetTitle( "r (cm)" );
        h->GetYaxis()->SetTitle( label[ih] );

        // rescale as rdphi
        if( use_phi_as_radian && ih ==1 )
        {
          for( int iy = 0; iy < h->GetNbinsY(); ++iy )
          {
            const double r = h->GetYaxis()->GetBinCenter( iy+1 );
            for( int ix = 0; ix < h->GetNbinsX(); ++ix )
            { h->SetBinContent( ix+1, iy+1, h->GetBinContent( ix+1, iy+1 )*r ); }
          }
        }

        cv->cd(ih+3*isuffix+1);
        gPad->SetRightMargin( 0.22 );

        h->Draw("hist");
        Draw::HorizontalLine( gPad, 30 )->Draw();
      }
    }

    pdfDocument.Add( cv );
  }

  return pdfFile;

}
