#include <RootUtil/Draw.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH3.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

static constexpr double z_rec = 5;


//_______________________________________________
TString DirectLaserEntries()
{
  set_style( false );
  gStyle->SetMarkerSize(1.);

  const TString tag = "_directlasers-new";
  const TString inputFile = "DST/TpcDirectLaserReconstruction.root";
  const auto pdfFile = Form( "Figures/DirectLaserEntries%s.pdf", tag.Data() );
  
  FileManager f( inputFile );
  PdfDocument pdfDocument( pdfFile );
  
  auto hentries = static_cast<TH3*>( f.GetHistogram( "entries" ) );
 
  // find matching zbin
  const auto zbin_rec = hentries->GetZaxis()->FindBin( z_rec );
  std::cout << "DirectLaserEntries - zbin_rec: " << zbin_rec << std::endl;
  
  {
    // entries in r, phi plane
    TCanvas* cv( new TCanvas( "cv1", "cv1", 800, 800 ) );
    cv->SetRightMargin( 0.24 );

    hentries->GetZaxis()->SetRange( zbin_rec, zbin_rec ); // z axis
    auto proj = hentries->Project3D( "yx" );
    proj->SetTitle("");
    proj->GetZaxis()->SetTitle( "entries" );
    proj->GetZaxis()->SetTitleOffset( 2.1 );
    proj->Draw( "colz" );  

    pdfDocument.Add(cv);

  }

  return pdfFile;
  
}
