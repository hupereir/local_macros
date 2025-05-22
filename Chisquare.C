#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

#include <TFile.h>
#include <TStyle.h>

#include "LayerDefines.h"

R__LOAD_LIBRARY(libRootUtilBase.so)

//__________________________________________
TString Chisquare()
{
  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

//   const TString tag = "_realistic_micromegas-new";
//   const TString inputFile = "DST/CONDOR_realistic_micromegas/dst_reco_truth_notpc_distortions-new/dst_reco_realistic_micromegas_*.root";

  const TString tag = "_realistic_micromegas";
  const TString inputFile = "DST/CONDOR_realistic_micromegas/dst_reco_truth_notpc_distortions/dst_reco_realistic_micromegas_*.root";

  // pdf output
  TString pdfFile( Form( "Figures/chisquare%s.pdf", tag.Data()  ) );
  PdfDocument pdfDocument( pdfFile );

  // get tree
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return pdfFile;
  
  const TCut cut = 
    "_tracks._nclusters_mvtx==3 "
    "&& _tracks._nclusters_intt>=2"
    "&& _tracks._nclusters_micromegas>=2";
  
  const TCut ndfCut = "_tracks._ndf >= 9";
  // const TCut cut;

  {
    const TString var = "_tracks._chisquare/_tracks._ndf";
    
    const TString cvName( "cv" );
    auto cv = new TCanvas( cvName, cvName, 800, 800 );
    auto legend = new TLegend( 0.5, 0.8, 0.97, 0.9, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    
    const TString hname( "chisquare" );
    // auto h = new TH1F( hname, "", 100, 0, 10 );
    auto h = new TH1F( hname, "", 100, 0, 1 );
    Utils::TreeToHisto( tree, hname, var, cut && ndfCut, false );
    h->SetTitle( "");
    h->GetXaxis()->SetTitle( "#chi^{2}" );
    
    h->Draw();
    gPad->SetLogy( true );
    
    legend->AddEntry( h, Form( "mean=%.3f#pm%.3f", h->GetMean(), h->GetMeanError() ), "PL" );
    
    std::cout << "Chisquare - Mean: " << h->GetMean() << " +/- " << h->GetMeanError() << std::endl;
    legend->Draw();
    
    pdfDocument.Add(cv);
  }
  
  {
    const TString var = "_tracks._ndf";
    
    const TString cvName( "cv" );
    auto cv = new TCanvas( cvName, cvName, 800, 800 );
    auto legend = new TLegend( 0.5, 0.8, 0.97, 0.9, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    
    const TString hname( "ndf" );
    auto h = new TH1F( hname, "", 20, 0, 20 );
    Utils::TreeToHisto( tree, hname, var, cut, false );
    h->SetTitle( "");
    h->GetXaxis()->SetTitle( "ndf" );
    
    h->Draw();
    
    legend->AddEntry( h, Form( "mean=%.3f#pm%.3f", h->GetMean(), h->GetMeanError() ), "PL" );
    
    std::cout << "Chisquare - Mean: " << h->GetMean() << " +/- " << h->GetMeanError() << std::endl;
    legend->Draw();
    
    pdfDocument.Add(cv);
  }
  
  return pdfFile;
}
