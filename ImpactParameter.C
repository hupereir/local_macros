#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
void ImpactParameter()
{
  
  set_style( false );
  
//   const TString tag = "_Hijing_Micromegas_merged" ;
//   const TString files[] = 
//   {
//     "DST/CONDOR_Hijing_Micromegas/dst_eval_merged/dst_eval_sHijing_0-12fm_merged_0*.root",
//     TString()
//   };

//   const TString tag = "_Hijing_Micromegas" ;
//   const TString files[] = 
//   {
//     "DST/CONDOR_Hijing_Micromegas/dst_eval/dst_eval_sHijing_0-12fm_0*.root",
//     TString()
//   };

  const TString tag = "_Hijing_Christof";
  const TString files[] = 
  {
    "DST/CONDOR_Hijing_Christof/dst_eval*.root",
    TString()
  };

  const TString pdfFile = Form( "Figures/ImpactParameter%s.pdf", tag.Data() );
  std::cout << "ClusterMultiplicity - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );
  
  // file manager
  FileManager fileManager;
  for( int i = 0; !files[i].IsNull(); ++i ) 
  { fileManager.AddFiles( files[i] ); }
  
  auto tree = fileManager.GetChain( "T" );
  
  if( true )
  {
    const TString var = "_events._bimp";
    auto h = new TH1F( "impact parameter", "", 100, 0, 13 );
    Utils::TreeToHisto( tree, h->GetName(), var, TCut(), false );
    h->GetXaxis()->SetTitle( "impact parameter b (fm)" );
    h->GetYaxis()->SetTitle( "dN/db" );
    
    // to get a flat distribution you need to normalize by the surface of the corresponding circle
    if( false )
    {
      for( int i = 0; i < h->GetNbinsX(); ++i )
      {
        auto b = h->GetXaxis()->GetBinCenter(i+1);
        auto entries = h->GetBinContent( i+1 );
        h->SetBinContent( i+1, entries/b );
      }
      h->GetYaxis()->SetTitle( "dN/db^{2}" );
    }
    
    std::unique_ptr<TCanvas> cv( new TCanvas( "cv", "cv", 800, 800 ) );
    h->Draw();
    pdfDocument.Add( cv.get() );
  }
  
  if( true )
  {
    const TString var = "100.*_events._bimp*_events._bimp/(12*12)";
    auto h = new TH1F( "centrality", "", 110, 0, 110 );
    Utils::TreeToHisto( tree, h->GetName(), var, TCut(), false );
    h->GetXaxis()->SetTitle( "centrality (%)" );
        
    std::unique_ptr<TCanvas> cv( new TCanvas( "cv1", "cv1", 800, 800 ) );
    h->Draw();
    pdfDocument.Add( cv.get() );
  }
  
}
