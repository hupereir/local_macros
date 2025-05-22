#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
// centrality from impact parameter
float get_centrality( float b )
{
  static constexpr float bmax=12;
  return 100.*::square(b/bmax); 
}

//____________________________________________________________________________
TString ClusterSize_Micromegas()
{
  
  set_style( false );
  const TString tag = "_realistic_full_micromegas_nominal";
  const TString files[] =
  {
    Form( "DST/CONDOR%s/dst_eval%s*.root", tag.Data(), tag.Data() ),
    TString()
  };

//   const TString tag = "_Hijing_Micromegas_merged" ;
//   const TString files[] = 
//   {
//     "DST/CONDOR_Hijing_Micromegas/dst_eval_merged-new/dst_eval_sHijing_0-12fm_merged_0*.root",
//     TString()
//   };

//   const TString tag = "_Hijing_Micromegas_single";
//   const TString files[] = 
//   {
//     "DST/CONDOR_Hijing_Micromegas/dst_eval/dst_eval_sHijing_0-12fm_0*.root",
//     TString()
//   };

  const TString rootFile = Form( "Rootfiles/ClusterSize_Micromegas%s.root", tag.Data() );
  const TString pdfFile = Form( "Figures/ClusterSize_Micromegas%s.pdf", tag.Data() );
  std::cout << "ClusterSize_Micromegas - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );
  
  // file manager
  FileManager fileManager;
  for( int i = 0; !files[i].IsNull(); ++i ) 
  { fileManager.AddFiles( files[i] ); }

  auto tree = fileManager.GetChain( "T" );

  // store histograms
  std::vector<TH1*> histograms;
  
  if( true )
  { 
    // cluster size per event per tile vs track angle
    // TCut cut( "get_centrality(DST#EVAL#SimEvaluator_hp::Container._events[0]._bimp)>70" );
    TCut cut;
    std::unique_ptr<TCanvas> cv( new TCanvas( "cv1", "cv1", 800, 400 ) );
    cv->Divide( 2, 1 );

    {
      TString pname = "profile_55";
      auto p = new TProfile( pname, pname, 100, -1, 1 );
      auto var = "_clusters[]._size:_clusters[]._truth_alpha";
      TCut layer_cut = Form( "_clusters[]._layer == %i", 55 );
      Utils::TreeToHisto( tree, pname, var, layer_cut&&cut, false );
      
      p->SetTitle("");
      p->GetXaxis()->SetTitle( "#alpha_{truth}" );
      p->GetYaxis()->SetTitle( Form( "csize (layer %i)", 55 ));
      
      cv->cd(1);
      p->SetMarkerStyle(20);
      p->SetMinimum(0);     
      // p->SetMaximum(5);
      p->Draw();
      
      histograms.push_back( p );
    }

    {
      TString pname = "profile_56";
      auto p = new TProfile( pname, pname, 100, -1, 1 );
      auto var = "_clusters[]._size:_clusters[]._truth_beta";
      TCut layer_cut = Form( "_clusters[]._layer == %i", 56 );
      Utils::TreeToHisto( tree, pname, var, layer_cut&&cut, false );
      
      p->SetTitle("");
      p->GetXaxis()->SetTitle( "#beta_{truth}" );
      p->GetYaxis()->SetTitle( Form( "csize (layer %i)", 56 ));
      
      cv->cd(2);
      p->SetMarkerStyle(20);
      p->SetMinimum(0);     
      // p->SetMaximum(5);
      p->Draw();

      histograms.push_back( p );
    }
    
    pdfDocument.Add( cv.get() );

  }
  
  if( false )
  {
    // cluster size per event per tile vs centrality
    std::unique_ptr<TCanvas> cv( new TCanvas( "cv1", "cv1", 800, 400 ) );
    cv->Divide( 2, 1 );
    
    for( int i = 0; i < 2; ++i )
    {
      // get multiplicity vs centrality
      TString pname = Form( "profile_%i", i );
      auto p = new TProfile( pname, pname, 10, 0, 100 );
      auto var = "DST#EVAL#TrackingEvaluator_hp::Container._clusters[]._size:get_centrality(DST#EVAL#SimEvaluator_hp::Container._events[0]._bimp)";
      TCut layer_cut = Form( "DST#EVAL#TrackingEvaluator_hp::Container._clusters[]._layer == %i", 55+i );
      TCut cut;
      // TCut cut( "DST#EVAL#TrackingEvaluator_hp::Container._clusters[]._truth_size==1" );
      Utils::TreeToHisto( tree, pname, var, layer_cut&&cut, false );

      p->SetTitle("");
      p->GetXaxis()->SetTitle( "Centrality (%)" );
      p->GetYaxis()->SetTitle( Form( "csize (layer %i)", 55+i ));
      
      cv->cd(i+1);
      p->SetMarkerStyle(20);
      p->SetMinimum(0);
      p->Draw();
      histograms.push_back( p );
      
    }
    
    pdfDocument.Add( cv.get() );
    
  }

  // save everything in rootfiles
  std::unique_ptr<TFile> output( TFile::Open( rootFile, "RECREATE" ) );
  output->cd();
  for( auto&& h:histograms) { if(h) h->Write(); }
  output->Close();

  
  return pdfFile;
      
}
