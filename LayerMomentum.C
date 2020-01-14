#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
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
void LayerMomentum( TString tag = TString() )
{

  const int layer = 54;
  const auto maxMomentum = 5;

  // pdf output
  if( tag.IsNull() ) tag = "_1k_truth_notpc_nominal" ;
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );
  const TString pdfFile = Form( "Figures/Momentum_truth%s_%i.pdf", tag.Data(), layer );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return;

  // variable names
  const TString var( "_clusters._truth_pt" );
  const TCut cut = Form( "_clusters._layer==%i && _clusters._truth_pt < %d", layer, maxMomentum );

  auto h = new TH1F( "h", "", 100, 0, maxMomentum );
  Utils::TreeToHisto( tree, h->GetName(), var, cut, kFALSE );

  auto cv = new TCanvas( "cv", "cv", 800, 800 );
  h->SetTitle("");
  h->GetXaxis()->SetTitle( "#it{p}_{T} (GeV/#it{c})" );
  h->Draw();

  Draw::PutText( 0.2, 0.8, Form( "layer = %i", layer ) );
  cv->SaveAs( pdfFile );
}
