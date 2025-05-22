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

//____________________________________________________________________________
float delta_phi( float phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
TString DeltaRPhi_QA_2d()
{

  set_style( false );

  // input files
  // const TString tag = "_TrackReconstruction";
  const TString tag = "_TrackReconstruction_genfit";
  const TString inputFile = Form( "Rootfiles/QA%s/DistortionsQA-00053877-*-full.root", tag.Data() );
  const TString pdfFile = Form( "Figures/DeltaRPhi_QA_2d%s.pdf", tag.Data() );

  std::cout << "DeltaRPhi - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhi - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto h = static_cast<TH2*>(fileManager.GetHistogram( "h_QAG4SimulationDistortions_deltarphi_layer" ));
  auto p = h->ProfileX();

  auto cv( new TCanvas( "cv", "cv", 900, 600 ) );
  cv->SetRightMargin(0.2);

  // h->SetStats(0);
  h->Draw("colz");
  p->SetLineColor(1);
  p->SetMarkerColor(1);
  p->SetMarkerStyle(20);
  p->Draw("same");

  gPad->SetLogz(true);
  gPad->Update();
  Draw::HorizontalLine(gPad, 0)->Draw();

  pdfDocument.Add( cv );

  return pdfFile;
}
