#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
void Radiograph()
{

  set_style( false );

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // pdf output
  const TString tag = "_realistic_micromegas_config3";
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval%s*.root", tag.Data(), tag.Data() );
  const TString pdfFile = Form( "Figures/Radiograph%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // configuration
  const Bool_t doFit = kTRUE;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // variable names
  {
    const TString var( "_clusters._truth_y:_clusters._truth_x" );
    const TString hName = "radiograph";
    auto h2d = new TH2F( hName, "", 500, -90, 90, 500, -90, 90 );
    Utils::TreeToHisto( tree, hName, var, TCut(), kFALSE );

    auto cv = new TCanvas( "", "", 800, 800 );
    h2d->GetXaxis()->SetTitle( "x (cm)" );
    h2d->GetYaxis()->SetTitle( "y (cm)" );
    h2d->SetTitle( "" );
    h2d->Draw( "col" );

    pdfDocument.Add( cv );
  }

  // variable names
  {
    const TString var( "_clusters._truth_r:_clusters._truth_z" );
    const TString hName = "radiograph2";
    auto h2d = new TH2F( hName, "", 500, -110, 110, 500, 0, 90 );
    Utils::TreeToHisto( tree, hName, var, TCut(), kFALSE );

    auto cv = new TCanvas( "", "", 800, 800 );
    h2d->GetXaxis()->SetTitle( "z (cm)" );
    h2d->GetYaxis()->SetTitle( "r (cm)" );
    h2d->SetTitle( "" );
    h2d->Draw( "col" );

    pdfDocument.Add( cv );
  }

}
