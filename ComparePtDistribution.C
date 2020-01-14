#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>

#include <TCanvas.h>
#include <TH1.h>

#include <array>

R__LOAD_LIBRARY(libRootUtilBase.so)

template< class T > T square( const T& x ) { return x*x; }

//___________________________________________
void ComparePtDistribution()
{

  set_style( false );

  constexpr std::array<float, 34> pt = {{0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.1, 2.3, 2.5, 2.7, 2.9, 3.12, 3.38, 3.62, 3.88, 4.2, 4.6, 5, 5.4, 5.8, 6.25, 6.75, 7.5, 8.5, 9.5}};
  constexpr std::array<float, 36> ptbins = {{0, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.5, 3.8, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.5, 7, 8, 9, 10}};
  constexpr std::array<float, 34> yield = {{22.3, 14.6, 9.76, 6.63, 4.57, 3.21, 2.29, 1.65, 1.19, 0.866, 0.628, 0.458, 0.337, 0.248, 0.183, 0.115, 0.0642, 0.0362, 0.0206, 0.0119, 0.00661, 0.00353, 0.00195, 0.0011, 0.000544, 0.000241, 0.000112, 6.07e-05, 3.06e-05, 1.58e-05, 8.85e-06, 4.05e-06, 1.45e-06, 9.38e-07}};
  constexpr std::array<float, 34> stat = {{0.0067, 0.0044, 0.0039, 0.0027, 0.0023, 0.0019, 0.0014, 0.0012, 0.00095, 0.00078, 0.00063, 0.00055, 0.00044, 0.00037, 0.00031, 0.00016, 0.00012, 8.3e-05, 6.2e-05, 4.5e-05, 2.9e-05, 2e-05, 1.5e-05, 1.1e-05, 5.8e-06, 3.8e-06, 2.5e-06, 1.9e-06, 1.4e-06, 1e-06, 8e-07, 4.4e-07, 3.9e-07, 3.6e-07}};
  constexpr std::array<float, 34> syst = {{1.7, 1.1, 0.72, 0.49, 0.34, 0.24, 0.17, 0.12, 0.088, 0.064, 0.046, 0.034, 0.025, 0.018, 0.014, 0.0085, 0.0048, 0.0027, 0.0015, 0.00088, 0.00052, 0.00028, 0.00015, 8.7e-05, 4.3e-05, 1.9e-05,8.8e-06, 4.8e-06, 2.7e-06, 1.6e-06, 1e-06, 5.2e-07, 3.2e-07, 3e-07}};
  constexpr std::array<float, 34> global = {{0.8, 0.53, 0.35, 0.24, 0.16, 0.12, 0.082, 0.059, 0.043, 0.031, 0.023, 0.016, 0.012, 0.0089, 0.0066, 0.0041, 0.0023, 0.0013, 0.00074, 0.00043, 0.00024, 0.00013, 7e-05, 4e-05, 2e-05, 8.7e-06, 4e-06, 2.2e-06, 1.1e-06, 5.7e-07, 3.2e-07, 1.5e-07, 5.2e-08, 3.4e-08}};

//   constexpr std::array<float,39> ptbins = {{0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.5, 3.8, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.5, 7, 8, 9, 10}};
//   constexpr std::array<float,37> yield = {{17.7, 10, 6.01, 3.6, 2.22, 1.44, 0.897, 0.596, 0.462, 0.33, 0.238, 0.171, 0.125, 0.0905, 0.066, 0.0485, 0.0357, 0.0264, 0.0166, 0.00925, 0.00521, 0.00297, 0.00171, 0.000952, 0.000509, 0.000281, 0.000158, 7.84e-05, 3.47e-05, 1.61e-05, 8.74e-06, 4.41e-06, 2.28e-06, 1.27e-06, 5.83e-07, 2.09e-07, 1.35e-07}};
//   constexpr std::array<float,37> error = {{0.0212, 0.00686, 0.00264, 0.00108, 0.000455, 0.000242, 0.000107, 5.43e-05, 0.0012, 0.0006, 0.000299, 0.000161, 8.5e-05, 4.39e-05, 2.4e-05, 1.3e-05, 6.73e-06, 4.07e-06, 1.5e-06, 4.78e-07, 1.51e-07, 4.68e-08, 1.61e-08, 5.63e-09, 1.64e-09, 4.72e-10, 1.6e-10, 3.91e-11, 7.79e-12, 1.74e-12, 5.53e-13, 1.92e-13, 7.39e-14, 3.4e-14, 9.63e-15, 5.28e-15, 4.56e-15}};

  const TString pdfFile( "Figures/ComparePtDistribution-new.pdf" );
  PdfDocument pdfDocument( pdfFile );

  // build histogram
  auto h = new TH1F( "h", "", ptbins.size()-1, &ptbins[0] );
  std::cout << "max: " << h->GetXaxis()->GetXmax() << std::endl;
  for( int i = 0; i < yield.size(); ++i )
  {
    h->SetBinContent( i+2, yield[i] );
    // h->SetBinError( i+2, error[i] );
    h->SetBinError( i+2, std::sqrt( square( stat[i] ) + square( syst[i] ) ) );
  }

  // open DST
  // const TString inputFile = "DST/dst_eval_1k_realistic2_truth_notpc_*.root";
  const TString inputFile = "DST/CONDOR_realistic_*/dst_eval_realistic_truth_notpc_*.root";
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return;

  // variable names
  const TString var( "_mc_tracks._pt" );
  auto h1 = new TH1F( "h1", "", 100, 0, 10 );
  Utils::TreeToHisto( tree, h1->GetName(), var, TCut(), false );

  auto cv = new TCanvas( "cv", "cv", 800, 800 );
  h->SetMarkerStyle(20);
  h->SetMarkerColor(1);
  h->Scale( 1.0/h->Integral() );
  h->Draw( "E" );
  h->GetXaxis()->SetTitle( "#it{p}_{T} (Gev/#it{c})" );

  h1->SetMarkerStyle(20);
  h1->SetMarkerColor(2);
  h1->SetLineColor(2);
  h1->Scale( 1.0/h1->Integral() );
  h1->Draw( "E Same" );


  auto legend = new TLegend( 0.29, 0.82, 0.94, 0.9, "", "NDC" );
  // legend->SetFillColor(0);
  // legend->SetFillStyle(0);
  // legend->SetBorderSize(0);
  legend->Draw();
  legend->AddEntry( h, "charged hadrons (PRC69 (2004) 034910)", "AP" );
  legend->AddEntry( h1, "#pi^{+} (this simulation)", "AP" );

  gPad->SetLogy( true );

  pdfDocument.Add( cv );

  // mean values
  std::cout << "ComparePtDistribution - mean " << h->GetMean() << "+/-" << h->GetMeanError() << " vs " << h1->GetMean() << "+/-" << h1->GetMeanError() << std::endl;

}
