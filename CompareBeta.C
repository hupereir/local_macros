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
#include <TStyle.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
TString CompareBeta()
{

  set_style( false );

  // initial guess for max residuals
  std::array<float, nDetectors> max_det_residual = { 0.003, 0.01, 0.2, 0.2, 0.2, 0.5, 0.5};

  // input files
  const TString tag = "_acts_truth_notpc_nodistortion-test";
  const TString inputFile = Form( "DST/CONDOR_realistic_micromegas/dst_reco%s/dst_reco_realistic_micromegas_1?.root", tag.Data() );


  const TString pdfFile = Form( "Figures/CompareBeta%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/CompareBeta%s.root", tag.Data() );

  std::cout << "CompareBeta - inputFile: " << inputFile << std::endl;
  std::cout << "CompareBeta - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool do_fit = true;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString prefix = "DST#EVAL#PHTpcResiduals::Container";
//   const TString var =
//     Form( "%s._tracks._clusters[]._trkr_acts._beta - %s._tracks._clusters[]._trkr_genfit._beta", 
//     prefix.Data(), prefix.Data() );
  const TString var =
    Form( "%s._tracks._clusters[]._trkr_acts._beta", 
    prefix.Data() );
  const TString var2d = Form( "%s:%s._tracks._clusters[]._trkr_genfit._beta", var.Data(), prefix.Data() );
  
  const double maxrphi = 0.1;
  const TCut cut;
  
  const double maxbeta= 1.;
  
  const TString hname = "deltabeta"; 
  TH2* h = new TH2F( hname, hname, 500, -maxbeta, maxbeta, 500, -maxbeta, maxbeta );
  Utils::TreeToHisto( tree, hname, var2d, cut, false );
  
  const TString cvName = "cv";
  auto cv( new TCanvas( cvName, cvName, 800, 800 ) );

  h->Draw();
  pdfDocument.Add( cv );
  
  return pdfFile;

}
