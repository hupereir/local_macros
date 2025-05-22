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


TString truth_pt()
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
  const TString pdfFile = Form( "Figures/truth_pt_all%s.pdf", global_tag.c_str());
  PdfDocument pdfDocument( pdfFile );

  // create canvas
  auto cv( new TCanvas( "cv", "cv", 1200, 800 ) );
  cv->Divide(4,2);

  int cvid = 0;
  for( const auto& [pid, tag]:data_list )
  {

    const TString inputFile = Form( "DST/CONDOR%s/DST_RECO%s/dst_reco*.root", tag.c_str(), global_tag.c_str() );

    FileManager fileManager( inputFile );
    auto tree = fileManager.GetChain( "T" );
    if( !tree ) { continue; }

    const TString var( "_tracks._truth_pt" );
    const TCut pid_cut = Form("_tracks._pid == %i", pid);

    auto h = Utils::TreeToHisto( tree, Form( "pt_%s", tag.c_str()), var, pid_cut, true );
    h->SetTitle(tag.c_str());
    h->SetFillStyle(1001);
    h->SetFillColor(kYellow );
    // h->SetStats(false);
    h->GetXaxis()->SetTitle( "#it{p}_{T,truth} (GeV)" );

    cv->cd(++cvid);
    h->Draw();
  }

  pdfDocument.Add(cv);

  return pdfFile;
}
