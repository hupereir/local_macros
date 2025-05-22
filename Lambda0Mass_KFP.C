#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

void Lambda0Mass_KFP()
{

  set_style(kFALSE);
  TGaxis::SetMaxDigits(3);

  const std::vector<int> runnumbers = {53877};
  // const std::vector<int> runnumbers = {53534,53631,53687,53716,53738,53739,53741,53742,53743,53744,53756,53783,53876,53877,53879};
  const auto runLabel = Utils::GetRunLabel(runnumbers);
  const auto runPostfix = Utils::GetRunPostfix(runnumbers);

  // load files
  FileManager filemanager;
  for( const auto& runnumber:runnumbers )
  {
    const auto inputFile = Form( "Rootfiles/KFP/lambda0/KFP-%08i-merged.root", runnumber );
    filemanager.AddFiles(inputFile);
  }

  auto tree = filemanager.GetChain( "DecayTree" );

  std::cout << "entries: " << tree->GetEntries() << std::endl;

  auto h = new TH1F( "lambda0mass", "", 110, 1.08, 1.19);
  double bin_width =(h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin())/h->GetNbinsX();
  h->GetXaxis()->SetTitle( "m(p^{#pm}#pi^{#mp}) [GeV]" );
  h->GetYaxis()->SetTitle( Form("candidates/(%.0f MeV)", bin_width*1000));

  // cuts
  const TCut mvtx_states_cut("min(track_1_MVTX_nStates,track_2_MVTX_nStates)>=2");
  const TCut tpc_states_cut("min(track_1_TPC_nStates,track_2_TPC_nStates)>25");
  const TCut track_chi2_cut(
    "(track_1_nDoF>0&&track_1_chi2/track_1_nDoF>=0&&track_1_chi2/track_1_nDoF<100)&&"
    "(track_2_nDoF>0&&track_2_chi2/track_2_nDoF>=0&&track_2_chi2/track_2_nDoF<100)"
    );
  const TCut track_pt_cut( "min(track_1_pT,track_2_pT)>=0.2");
  const TCut track_ip_cut( "min(fabs(track_1_IP_xy), fabs(track_2_IP_xy))>0.05");
  const TCut track_dca_cut("track_1_track_2_DCA<0.5");

  const TCut cut= mvtx_states_cut && tpc_states_cut&&track_chi2_cut
    &&track_pt_cut
    &&track_ip_cut
    &&track_dca_cut
    &&TCut();

  Utils::TreeToHisto( tree, h->GetName(), "Lambda0_mass", cut, false );

  // create canvas
  auto cv( new TCanvas( "cv", "cv", 800, 800 ) );
  h->Draw();
}
