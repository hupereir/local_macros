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

#include "LayerDefines.h"
#include "Fit.C"

R__LOAD_LIBRARY(libRootUtilBase.so)

TString VertexTime()
{

  const TString tag = "_realistic_micromegas_50khz";
  const TString inputFile = "DST/dst_g4hits_merged.root";

//   const TString tag = "_hijing_micromegas_50khz";
//   const TString inputFile = "DST/CONDOR_Hijing_Micromegas_50kHz/dst_simeval_merged/3/dst*.root";

  const TString pdfFile = Form( "Figures/VertexTime%s.pdf", tag.Data() );

  std::cout << "VertexTime - inputFile: " << inputFile << std::endl;
  std::cout << "VertexTime - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  {
    const TString var( "_vertex_list._t" );
    const TCut cut( "_vertex_list._is_main_vertex==0" );
    auto h = Utils::TreeToHisto( tree, "h_t", var, cut, true );
    h->SetTitle( "" );
    h->GetXaxis()->SetTitle( "#Delta t (ns)" );

    auto cv = new TCanvas( "cv", "cv", 900, 600 );
    h->Draw();
    pdfDocument.Add( cv );

  }


  {
    const TString var( "_events._nevt_bg" );
    const TCut cut;
    auto h = new TH1F( "h_t", "h_t", 10, -0.5, 9.5 );
    Utils::TreeToHisto( tree, "h_t", var, cut, false );
    h->SetTitle( "" );
    h->GetXaxis()->SetTitle( "Background events" );

    auto cv = new TCanvas( "cv1", "cv", 900, 600 );
    h->Draw();
    pdfDocument.Add( cv );

    std::cout << "VertexTime - mean background : " << h->GetMean() << std::endl;

  }

  return pdfFile;
}
