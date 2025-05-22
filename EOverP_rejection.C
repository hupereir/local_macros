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

//_____________________________________________________________________________
TH1* rejection( TH1* source )
{
  auto h = static_cast<TH1*>(source->Clone(Form( "%s_rejection", source->GetName())));
  h->Reset();
  h->GetXaxis()->SetTitle(source->GetXaxis()->GetTitle());
  h->GetYaxis()->SetTitle("rejection");

  const auto integral = source->Integral();
  for( int ibin=0; ibin<source->GetNbinsX(); ++ibin )
  {
    const auto bin_integral = source->Integral(ibin+1,source->GetNbinsX()+1);
    const auto ratio = bin_integral/integral;
    const auto error = std::sqrt(ratio*(1.-ratio)/integral);
    const auto ratio_inv = 1./ratio;
    const auto error_inv = error/square(ratio);
    h->SetBinContent(ibin+1, ratio_inv);
    h->SetBinError(ibin+1, error_inv);
  }
  return h;
}

//_____________________________________________________________________________
TString EOverP_rejection()
{

  set_style( false );
  // gStyle->SetOptStat(0);

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
  const std::string subtag = "_p2";

  // pdf document
  const TString pdfFile = Form( "Figures/e_over_p_rejection%s%s.pdf", global_tag.c_str(),subtag.c_str());
  PdfDocument pdfDocument( pdfFile );

  // rootfile
  const TString rootFilename = Form( "Rootfiles/e_over_p_rejection%s%s.root", global_tag.c_str(),subtag.c_str());
  RootFile rootFile( rootFilename );

  // Utils::max_entries = 10000;

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

    const TString var( "_tracks._calo_clusters._e/_tracks._p" );
    const TCut cut = "_tracks._calo_clusters._layer == 1";
    const TCut pt_cut = "_tracks._p>2";
    const TCut pid_cut = Form("_tracks._pid == %i", pid);

    auto h = new TH1F( Form("e_over_p%s", tag.c_str()), "", 100, 0, 3 );
    Utils::TreeToHisto( tree, h->GetName(), var, cut&&pt_cut&&pid_cut, false );
    h->SetTitle(tag.c_str());
    h->GetXaxis()->SetTitle( "E/p" );

    auto h_rejection = rejection(h);
    h_rejection->SetMarkerStyle(20);
    h_rejection->SetMarkerColor(1);
    h_rejection->SetStats(0);

    cv->cd(++cvid);
    h_rejection->Draw();

    rootFile.Add(h);
    rootFile.Add(h_rejection);

  }

  pdfDocument.Add(cv);

  return pdfFile;
}
