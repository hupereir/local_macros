#include <RootUtil/PdfDocument.h>
#include <TString.h>
#include <TFile.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

TString DrawDistortions()
{
  set_style( false );
  gStyle->SetPadLeftMargin(0.13);

  // load tfile
  const TString file = "PHSpaceChargeReconstruction.root";
  std::unique_ptr<TFile> tfile( TFile::Open( file ) );

  const TString pdfFile = "Figures/PHSpaceChargeReconstruction.pdf";
  PdfDocument pdfDocument( pdfFile );

  // reconstructed
  static constexpr int ncoord = 3;
  std::array<TGraphErrors*, ncoord> tg;
  std::array<TGraphErrors*, ncoord> tg_fit;
  std::array<TString, ncoord> y_label = { "r#Delta#phi (cm)", "#Deltaz (cm)", "#Deltar (cm)" };
  for( int i = 0; i < ncoord; ++i )
  {

    // fit
    constexpr int iz = 0;
    constexpr int iphi = 0;
    TString tgname = Form( "tg_%i_%i_%i", iz, iphi, i );
    tg_fit[i] = static_cast<TGraphErrors*>( tfile->Get( tgname ) );

    if( !tg_fit[i] )
    {
      std::cout << "DrawDistortions - unable to find " << tgname << std::endl;
      continue;
    }

    // reference
    tg[i] = new TGraphErrors();


    // loop over tg_fit entries
    for( int ipoint=0; ipoint < tg_fit[i]->GetN(); ++ipoint )
    {
      double r, y;
      tg_fit[i]->GetPoint( ipoint, r, y );

      if( i==2 && true )
      {

        static constexpr float deltar_max = 1.5;
        const float dr = -deltar_max*std::cos( 2*M_PI*(r-rmin_tpc)/rlength_tpc );
        tg[i]->SetPoint( ipoint, r, dr );

      } else tg[i]->SetPoint( ipoint, r, 0 );

    }

    // make plot
    auto cv_name = Form( "cv%i", i );
    auto cv( new TCanvas( cv_name, cv_name, 800, 800 ) );
    tg[i]->SetMarkerStyle( 20 );
    tg[i]->SetMarkerColor( 4 );
    tg[i]->SetLineColor( 4 );
    tg[i]->GetXaxis()->SetTitle( "r (cm)" );
    tg[i]->GetYaxis()->SetTitle( y_label[i] );
    tg[i]->SetMinimum( -2.5 );
    tg[i]->SetMaximum( 2.5 );
    tg[i]->Draw("AP");

    tg_fit[i]->SetMarkerStyle( 20 );
    tg_fit[i]->SetMarkerColor( 2 );
    tg_fit[i]->SetLineColor( 2 );
    tg_fit[i]->Draw("P");
    cv->Update();

    auto legend = new TLegend( 0.16, 0.82, 0.60, 0.93, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();
    legend->AddEntry( tg[i], "input distortion", "AP" );
    legend->AddEntry( tg_fit[i], "reconstructed distortion", "AP" );

    pdfDocument.Add( cv );

  }

  return pdfFile;
}
