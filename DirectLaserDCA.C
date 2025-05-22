#include <RootUtil/Draw.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH3.h>

R__LOAD_LIBRARY(libRootUtilBase.so)


//_______________________________________________
TString DirectLaserDCA()
{
  set_style( false );
  gStyle->SetMarkerSize(1.);

  const TString tag = "_all_directlasers_simple_no_distortion-new";
  const TString inputFile = Form( "DST/TpcDirectLaserReconstruction%s.root", tag.Data() );
  const auto pdfFile = Form( "Figures/DirectLaserDCA%s.pdf", tag.Data() );
  
  FileManager f( inputFile );
  PdfDocument pdfDocument( pdfFile );
  
  auto hdca = static_cast<TH2*>( f.GetHistogram( "dca_layer" ) );
   
  {
    // entries in r, phi plane
    TCanvas* cv( new TCanvas( "cv1", "cv1", 800, 800 ) );

    auto proj = hdca->ProjectionY();
    proj->SetTitle("");
    proj->GetXaxis()->SetTitle( "DCA (cm)" );
    proj->GetZaxis()->SetTitleOffset( 2.1 );
    proj->Draw();  

    gPad->SetLogy( true );
    
    pdfDocument.Add(cv);

  }

  return pdfFile;
  
}
