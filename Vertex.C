#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TH1.h>

#include <array>

R__LOAD_LIBRARY(libRootUtilBase.so)

TString Vertex()
{

  set_style( false );
  const TString pdfFile( "Figures/Vertex_hijing_Micromegas.pdf" );
  PdfDocument pdfDocument( pdfFile );

  // open DST
  const TString inputFile = "DST/CONDOR_Hijing_Micromegas/dst_eval_*.root";
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();


  const TCut cut( "_vertex_list._is_main_vertex" );
  const TString var( "_vertex_list._z" );
  auto h = new TH1F( "h", "", 100, -100, 100 );
  Utils::TreeToHisto( tree, h->GetName(), var, cut, false );
  h->SetTitle( "" );

  auto cv = new TCanvas( "cv", "cv", 800, 800 );
  h->GetXaxis()->SetTitle( "z_{vtx} (cm)" );
  h->Draw();
  h->Fit( "gaus", "0" );
  auto f = h->GetFunction( "gaus" );
  f->SetLineColor(2);
  f->Draw( "same" );
  pdfDocument.Add( cv );

  return pdfFile;
}
