#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
TString ClusterMultiplicity_tpc()
{
  
  set_style( false );
    
  const TString tag = "_Hijing_Micromegas_merged" ;
  const TString files[] = 
  {
    "DST/CONDOR_Hijing_Micromegas/dst_eval_merged/dst_eval_sHijing_0-12fm_merged_0*.root",
    TString()
  };

//   const TString tag = "_Hijing_Micromegas_single";
//   const TString files[] = 
//   {
//     "DST/CONDOR_Hijing_Micromegas/dst_eval/dst_eval_sHijing_0-12fm_0*.root",
//     TString()
//   };
// 
//   const TString tag = "_Hijing_Christof";
//   const TString files[] = 
//   {
//     "DST/CONDOR_Hijing_Christof/dst_eval*.root",
//     TString()
//   };
  
  const TString pdfFile = Form( "Figures/ClusterMultiplicity%s_tpc.pdf", tag.Data() );
  std::cout << "ClusterMultiplicity_tpc - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager;
  for( int i = 0; !files[i].IsNull(); ++i ) 
  { fileManager.AddFiles( files[i] ); }
  
  auto tree = fileManager.GetChain( "T" );

  // get tpc cluster distribution
  const TString clus_var = "DST#EVAL#TrackingEvaluator_hp::Container._events._nclusters_tpc";
  auto hdist = Utils::TreeToHisto( tree, "nclusters_tpc", clus_var, TCut(), true );
  hdist->SetTitle( "" );
  hdist->GetXaxis()->SetTitle( "N_{clusters, TPC}" );
  
  {
    const TString cvName = "cv";
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
    cv->SetRightMargin( 0.1 );
    hdist->Draw();
    
    Draw::PutText( 0.4, 0.8, Form( "Mean = %.0f #pm %.0f", hdist->GetMean(), hdist->GetMeanError() ) );
    
    pdfDocument.Add( cv.get() );
  }

  return pdfFile;
  
}
