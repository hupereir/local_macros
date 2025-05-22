#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

void PileUp()
{

  set_style( false );

  const TString tag = "_Hijing_Micromegas_merged";
  const TString file = "DST/CONDOR_Hijing_Micromegas/dst_simeval_merged/dst_simeval_sHijing_0-12fm_merged_*.root";

//   const TString tag = "_Hijing_Micromegas_single";
//   const TString file = "DST/CONDOR_Hijing_Micromegas/dst_eval/dst_eval_sHijing_0-12fm_*.root";

//   const TString tag = "_Hijing_Christof";
//   const TString file = "DST/CONDOR_Hijing_Christof/dst_eval*.root";

  const TString pdfFile = Form( "Figures/PileUp%s.pdf", tag.Data() );
  std::cout << "PileUp - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager(file);
  auto tree = fileManager.GetChain( "T" );
  
  {
    auto h = Utils::TreeToHisto( tree, "h0", "_events._nevt_bg", TCut(), true );
    const TString cvName = "cv";
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    h->SetTitle( "" );
    h->Draw();
    Draw::PutText( 0.6, 0.8, Form( "Mean = %.2f #pm %.2f", h->GetMean(), h->GetMeanError() ) );
    pdfDocument.Add( cv.get() );
  }

  {
    // particle transverse momentum
    auto h = Utils::TreeToHisto( tree, "h0", "_particle_list[]._pt", "_particle_list[]._embed <= 0", true );

    const TString cvName = "cv1";
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    h->SetTitle( "" );
    h->Draw();
    cv->SetLogy( true );
    pdfDocument.Add( cv.get() );    
    
    std::cout << "PileUp - particles per event: " << h->GetEntries()/tree->GetEntries() << std::endl;
    
  }
  
}
