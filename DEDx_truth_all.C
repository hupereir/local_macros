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

#include <trackbase/TrkrClusterContainer.h>

R__LOAD_LIBRARY(libRootUtilBase.so)


TString dedx_truth_all()
{

  set_style( false );
  gStyle->SetOptStat(0);

  using data_type_pair_t = std::pair<int, std::string>;
  using data_list_t = std::vector<data_type_pair_t>;
  data_list_t data_list =
  {
    { 11, "_single_electron" },
    { -11, "_single_positron" },
    { 211, "_single_piplus" },
    { -211, "_single_piminus" },
    { 321, "_single_kplus" },
    { -321, "_single_kminus" },
    { 2212, "_single_proton" },
    { -2212, "_single_antiproton" }
  };

  const std::string global_tag = "_lpt";

  // pdf document
  const TString pdfFile = Form( "Figures/dedx_truth_all%s.pdf", global_tag.c_str());
  PdfDocument pdfDocument( pdfFile );

  // create dummy histogram
  auto h_all = new TH2F( "dedx", "", 100, -5, 5, 100, 0, 10 );
  h_all->GetXaxis()->SetTitle( "q.p_{truth} [GeV]" );
  h_all->GetYaxis()->SetTitle( "dE/dx_{truth} [eV/cm]" );
  h_all->GetYaxis()->SetTitleOffset(1.6);

  for( const auto& [pid, tag]:data_list )
  {

    const TString inputFile = Form( "DST/CONDOR%s/DST_RECO%s/dst_reco*.root", tag.c_str(), global_tag.c_str() );

    FileManager fileManager( inputFile );
    auto tree = fileManager.GetChain( "T" );
    if( !tree ) { continue; }

    const TString var( "_tracks._truth_dedx:_tracks._truth_p*_tracks._charge" );
    const TCut cut = "_tracks._truth_p < 5";
    const TCut pid_cut = Form("_tracks._pid == %i", pid);

    auto h = static_cast<TH2*>( h_all->Clone( Form("dedx_%s", tag.c_str() ) ) );
    h->Reset();
    Utils::TreeToHisto( tree, h->GetName(), var, cut && pid_cut, false );
    h_all->Add(h);
  }

  {
    // create canvas
    auto cv( new TCanvas( "cv", "cv", 800, 800 ) );
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.15);

    h_all->Draw("colz");
    gPad->SetLogz(true);
    pdfDocument.Add(cv);
  }

  {
    // create canvas
    auto cv( new TCanvas( "cv1", "cv1", 800, 800 ) );
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.15);

    h_all->Draw("colz");
    pdfDocument.Add(cv);
  }

  return pdfFile;
}
