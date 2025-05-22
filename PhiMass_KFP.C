#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>
#include <RootUtil/FitUtils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

static constexpr double mass_k = 0.493;
static constexpr double mass_phi = 1.019;

//____________________________________________________
double square(double x ) { return x*x; }

//____________________________________________________
double fit_function( double* x, double* par )
{
  const double mass = x[0];

  // gaussian signal
  const double amplitude = par[0];
  const double mean = par[1];
  const double sigma = par[2];
  const double signal = amplitude* UTILS::FitUtils::Gaus(mass, mean, sigma );

  const double mass_2k = par[3];
  if( mass < mass_2k )
  {
    TF1::RejectPoint();
    return 0;
  }

  const double amplitude_bg0 = par[4];
  const double amplitude_bg1 = par[5];
  const double background = std::sqrt(mass-mass_2k)*(amplitude_bg0 + amplitude_bg1*(mass-mass_2k));
  return signal+background;
}

//____________________________________________________
double fit_function_polynom( double* x, double* par )
{
  const double mass = x[0];

  // gaussian signal
  const double amplitude = par[0];
  const double mean = par[1];
  const double sigma = par[2];
  const double signal = amplitude* UTILS::FitUtils::Gaus(mass, mean, sigma );

  const double delta_mass = mass - 0.987;
  const double background = par[3] + par[4]*delta_mass + par[5]*square(delta_mass);
  return signal+background;
}

//____________________________________________________
void PhiMass_KFP()
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
    const auto inputFile = Form( "Rootfiles/KFP/phi/KFP-%08i-merged.root", runnumber );
    filemanager.AddFiles(inputFile);
  }

  auto tree = filemanager.GetChain( "DecayTree" );

  std::cout << "entries: " << tree->GetEntries() << std::endl;

  const double mass_min = 0.98;
  const double mass_max = 1.10;
  const int nbins = 120;
  const double bin_width = (mass_max-mass_min)/nbins;

  auto h = new TH1F( "phimass", "", nbins, mass_min, mass_max);
  h->GetXaxis()->SetTitle( "m(K^{+}K^{-}) [GeV]" );
  h->GetYaxis()->SetTitle( Form("candidates/(%.0f MeV)", bin_width*1000));

  // cuts
  const TCut mvtx_states_cut("min(track_1_MVTX_nStates,track_2_MVTX_nStates)>=2");
  const TCut tpc_states_cut("min(track_1_TPC_nStates,track_2_TPC_nStates)>25");
  const TCut track_chi2_cut(
    "(track_1_nDoF>0&&track_1_chi2/track_1_nDoF>=0&&track_1_chi2/track_1_nDoF<100)&&"
    "(track_2_nDoF>0&&track_2_chi2/track_2_nDoF>=0&&track_2_chi2/track_2_nDoF<100)"
    );
  const TCut track_pt_cut("min(track_1_pT,track_2_pT)>=0.2");
  const TCut track_p_cut("max(track_1_p,track_2_p)<=0.7");
  const TCut track_dca_cut("track_1_track_2_DCA<0.04");
  const TCut decay_length_cut("fabs(phi_decayLength)<0.05");

  // cut
  // const TCut cut = mvtx_states_cut&&track_pt_cut&&track_p_cut&&track_dca_cut&&decay_length_cut;
  const TCut cut= mvtx_states_cut && tpc_states_cut&&track_chi2_cut
    &&track_pt_cut
    &&track_p_cut
    &&track_dca_cut
    &&decay_length_cut;

  Utils::TreeToHisto( tree, h->GetName(), "phi_mass", cut, false );

  // create canvas
  auto cv( new TCanvas( "cv", "cv", 800, 800 ) );
  h->Draw();

//   // perform fit
//   auto f =  new TF1( "fit_function", fit_function, mass_min, mass_max, 6 );
//   f->SetParameter(0, h->GetMaximum());
//   f->SetParameter(1, mass_phi );
//   f->SetParameter(2, 0.002 );
//   f->SetParameter(3, 2*mass_k );
//   f->SetParameter(4, 0 );
//   f->SetParameter(5, 0 );

  // perform fit
  auto f =  new TF1( "fit_function", fit_function_polynom, mass_min, mass_max, 6 );
  f->SetParameter(0, h->GetMaximum());
  f->SetParameter(1, mass_phi );
  f->SetParameter(2, 0.002 );
  f->SetParameter(3, -15 );
  f->SetParameter(4, -0.4 );
  f->SetParameter(5, 15 );

  h->Fit( f, "0" );
  f->SetLineColor(1);
  f->Draw("same");

}
