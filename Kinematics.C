#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TH1.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
void Kinematics()
{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  set_style( false );

  // pdf output
  const TString tag = "_5k_truth_low_momentum" ;
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );
  const TString pdfFile = Form( "Figures/Kinematics%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  std::cout << "Entries: " << tree->GetEntries() << std::endl;

  {

    // create canvas
    const TString cvName = "cv_pt";
    auto cv = new TCanvas( cvName, cvName, 800, 800 );

    const TString var( "_mc_tracks._pt" );
    const TString hName( "pt" );
    std::unique_ptr<TH1> h( new TH1F( hName, "", 100, 0, 30 ) );
    Utils::TreeToHisto( tree, hName, var, TCut(), kFALSE );

    h->SetTitle( "" );
    h->Draw();
    h->GetXaxis()->SetTitle( "#it{p}_{T} (GeV/#it{c})" );
    pdfDocument.Add( cv );
  }

  {

    // create canvas
    const TString cvName = "cv_eta";
    auto cv = new TCanvas( cvName, cvName, 800, 800 );

    const TString var( "_mc_tracks._eta" );
    const TString hName( "eta" );
    std::unique_ptr<TH1> h( new TH1F( hName, "", 100, -1.5, 1.5 ) );
    Utils::TreeToHisto( tree, hName, var, TCut(), kFALSE );

    h->SetTitle( "" );
    h->Draw();
    h->GetXaxis()->SetTitle( "#eta" );
    pdfDocument.Add( cv );
  }

}
