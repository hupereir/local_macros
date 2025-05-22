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

TString DrawDistortionMap_cm()
{

  set_style( false );
  gStyle->SetOptStat(0);

  const TString tag = "_average";
  const TString inputfile = "/sphenix/tg/tg01/jets/bkimelman/BenProduction/Feb21_2025/Laminations_run2pp_ana464_2024p011_v001-00053534.root";
  const TString pdfFile( Form( "Figures/DistortionMap_cm%s.pdf", tag.Data() ) );
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
    "r#Delta#phi (cm)"
  };

  const TString cvname = "cv";
  TCanvas* cv( new TCanvas( cvname, cvname, 1200, 800 ) );
  cv->Divide( 3, 2 );

  for( int isuffix = 0; isuffix < 2; ++isuffix )
  {

    for( int ih = 0; ih < 3; ++ih )
    {
      const auto hname = Form( "%s%s", prefix[ih].Data(), suffix[isuffix].Data() );

      // deltaR vs R and z
      auto h= dynamic_cast<TH2*>(f->Get(hname));
      h->GetYaxis()->SetRangeUser( 30, h->GetYaxis()->GetXmax() );
      h->SetTitle( "" );

      h->GetXaxis()->SetTitle( "#phi (rad)" );
      h->GetYaxis()->SetTitle( "r (cm)" );
      h->GetZaxis()->SetTitle( label[ih] );
      h->GetZaxis()->SetTitleOffset( 1.8 );

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
      gPad->Update();
      Draw::HorizontalLine( gPad, 30 )->Draw();
    }
  }

  pdfDocument.Add( cv );
  return pdfFile;

}
