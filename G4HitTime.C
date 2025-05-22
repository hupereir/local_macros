#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
TString G4HitTime()
{

  static const float drift_velocity = 8./1000; // cm/ns
  static const float tpc_length = 211; // cm
  const auto tmax = (tpc_length/2)/drift_velocity;
  std::cout << "G4HitTime - tmax: " << tmax << std::endl;

  const TString tag = "_realistic_micromegas_50khz";
  const TString inputFile = "DST/dst_g4hits_merged.root";

  // const TString pdfFile = Form( "Figures/G4HitTime_tpc_single%s.pdf", tag.Data() );
  const TString pdfFile = Form( "Figures/G4HitTime_tpc%s.pdf", tag.Data() );

  std::cout << "G4HitTime - pdfFile: " << pdfFile << std::endl;
  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager ( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // TString var = Form( "(_g4hits._t)/1000", tpc_length/2, drift_velocity );
  TString var = Form( "(_g4hits._t+(%.2f-fabs(_g4hits._z))/%.4f)/1000", tpc_length/2, drift_velocity );
  TCut cut( "_g4hits._detid == 2&&_g4hits._embed<0" );
  // TCut cut( "_g4hits._detid == 2" );

  auto h = new TH1F( "h", "", 100, -1.1*tmax/1000, 2.*1.1*tmax/1000 );
  Utils::TreeToHisto( tree, h->GetName(), var, cut, false );
  h->SetTitle( "" );
  h->GetXaxis()->SetTitle("t_{0} + t_{drift} (#mus)");
  // h->GetXaxis()->SetMaxDigits(1);

  {
    const TString cvName = "cv";
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 600, 600 ) );
    cv->SetRightMargin( 0.1 );
    h->Draw();
    Draw::VerticalLine( cv.get(), 0 )->Draw();
    Draw::VerticalLine( cv.get(), tmax/1000 )->Draw();
    pdfDocument.Add( cv.get() );
  }

  // calculate integrals
  auto integral = h->Integral();
  const auto ibin_min = h->GetXaxis()->FindBin( 0. );
  const auto ibin_max = h->GetXaxis()->FindBin( tmax/1000 );
  auto accepted = h->Integral( ibin_min, ibin_max );
  std::cout << "G4HitTime - entries: " << h->GetEntries() << " integral: " << integral << " accepted: " << accepted << " ratio: " << accepted/integral << std::endl;

  return pdfFile;

}
