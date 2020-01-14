#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>

#include <TCanvas.h>
#include <TH1.h>

#include <array>

R__LOAD_LIBRARY(libRootUtilBase.so)

template< class T > T square( const T& x ) { return x*x; }

//___________________________________________
void ChargedHadronSpectra()
{

  set_style( false );

  constexpr std::array<float, 34> pt = {{0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.1, 2.3, 2.5, 2.7, 2.9, 3.12, 3.38, 3.62, 3.88, 4.2, 4.6, 5, 5.4, 5.8, 6.25, 6.75, 7.5, 8.5, 9.5}};
  constexpr std::array<float, 35> ptbins = {{0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.5, 3.8, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.5, 7, 8, 9, 10}};
  constexpr std::array<float, 34> yield = {{22.3, 14.6, 9.76, 6.63, 4.57, 3.21, 2.29, 1.65, 1.19, 0.866, 0.628, 0.458, 0.337, 0.248, 0.183, 0.115, 0.0642, 0.0362, 0.0206, 0.0119, 0.00661, 0.00353, 0.00195, 0.0011, 0.000544, 0.000241, 0.000112, 6.07e-05, 3.06e-05, 1.58e-05, 8.85e-06, 4.05e-06, 1.45e-06, 9.38e-07}};
  constexpr std::array<float, 34> stat = {{0.0067, 0.0044, 0.0039, 0.0027, 0.0023, 0.0019, 0.0014, 0.0012, 0.00095, 0.00078, 0.00063, 0.00055, 0.00044, 0.00037, 0.00031, 0.00016, 0.00012, 8.3e-05, 6.2e-05, 4.5e-05, 2.9e-05, 2e-05, 1.5e-05, 1.1e-05, 5.8e-06, 3.8e-06, 2.5e-06, 1.9e-06, 1.4e-06, 1e-06, 8e-07, 4.4e-07, 3.9e-07, 3.6e-07}};
  constexpr std::array<float, 34> syst = {{1.7, 1.1, 0.72, 0.49, 0.34, 0.24, 0.17, 0.12, 0.088, 0.064, 0.046, 0.034, 0.025, 0.018, 0.014, 0.0085, 0.0048, 0.0027, 0.0015, 0.00088, 0.00052, 0.00028, 0.00015, 8.7e-05, 4.3e-05, 1.9e-05,8.8e-06, 4.8e-06, 2.7e-06, 1.6e-06, 1e-06, 5.2e-07, 3.2e-07, 3e-07}};
  constexpr std::array<float, 34> global = {{0.8, 0.53, 0.35, 0.24, 0.16, 0.12, 0.082, 0.059, 0.043, 0.031, 0.023, 0.016, 0.012, 0.0089, 0.0066, 0.0041, 0.0023, 0.0013, 0.00074, 0.00043, 0.00024, 0.00013, 7e-05, 4e-05, 2e-05, 8.7e-06, 4e-06, 2.2e-06, 1.1e-06, 5.7e-07, 3.2e-07, 1.5e-07, 5.2e-08, 3.4e-08}};

  constexpr std::array<float,28> pt2 = {{0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95}};
  constexpr std::array<float,29> ptbins2 = {{0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3}};
  constexpr std::array<float,28> yield2 = {{107, 60.6, 36.3, 21.8, 13.4, 8.71, 5.41, 3.59, 2.35, 1.58, 1.05, 0.759, 0.516, 0.337, 0.244, 0.177, 0.127, 0.0901, 0.0668, 0.0471, 0.0327, 0.026, 0.0194, 0.0149, 0.0113, 0.0093, 0.0062, 0.00517}};
  constexpr std::array<float,28> stat2 = {{0.88, 0.5, 0.31, 0.2, 0.13, 0.095, 0.063, 0.045, 0.031, 0.022, 0.015, 0.012, 0.0083, 0.0056, 0.0042, 0.0033, 0.0024, 0.0019, 0.0012, 0.00089, 0.00068, 0.00062, 0.00053, 0.00047, 0.00042, 0.0004, 0.00032, 0.00031}};

  const TString pdfFile( "Figures/ChargedHadronSpectra.pdf" );
  PdfDocument pdfDocument( pdfFile );

  // normalization range
  double pt_min = 0.5;
  double pt_max = 3;

  // build first histogram
  std::vector<float> yield_integral;
  auto h = new TH1F( "h", "", ptbins.size()-1, &ptbins[0] );
  double_t norm = 0;
  for( int i = 0; i < yield.size(); ++i )
  {
    h->SetBinContent( i+1, yield[i] );
    h->SetBinError( i+1, std::sqrt( square( stat[i] ) + square( syst[i] ) ) );

    const auto binWidth = ptbins[i+1]-ptbins[i];
    yield_integral.push_back( yield[i]*binWidth );

    // normalization
    if( ptbins[i] >= pt_min && ptbins[i+1] <= pt_max ) norm += yield[i]*binWidth;

  }

  Stream::PrintVector<float>( "float", "yield_int", yield_integral, "%.3g" );
  h->Scale( 1.0/norm );

  // build second histogram
  std::vector<float> yield_integral2;
  auto h2 = new TH1F( "h2", "", ptbins2.size()-1, &ptbins2[0] );
  double norm2 = 0;
  for( int i = 0; i < yield2.size(); ++i )
  {
    h2->SetBinContent( i+1, yield2[i] );
    h2->SetBinError( i+1, stat2[i] );

    const auto binWidth = ptbins2[i+1]-ptbins2[i];
    yield_integral2.push_back( yield2[i]*binWidth );

    // normalization
    if( ptbins2[i] >= pt_min && ptbins2[i+1] <= pt_max ) norm2 += yield2[i]*binWidth;
  }

  Stream::PrintVector<float>( "float", "yield_int2", yield_integral2, "%.3g" );
  h2->Scale( 1.0/norm2 );

  // combine
  constexpr std::array<float,38> ptbins3 = {{
    0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1,
    1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2,
    2.4, 2.6, 2.8, 3, 3.2, 3.5, 3.8, 4, 4.4, 4.8,
    5.2, 5.6, 6, 6.5, 7, 8, 9, 10}};
  std::vector<float> yield3;
  std::vector<float> yield_integral3;
  std::vector<float> error3;

  auto h3 = new TH1F( "h3", "", ptbins3.size()-1, &ptbins3[0] );
  for( int i = 0; i < ptbins3.size()-1; ++i )
  {

    double y = 0;
    double sumw = 0;

    // add first histogram
    for( int i1 = 0; i1 < ptbins.size()-1; ++i1 )
    {
      if( ptbins3[i] >= ptbins[i1] && ptbins3[i+1] <= ptbins[i1+1] )
      {
        double w = 1.0/(square(stat[i1]/norm) + square(syst[i1]/norm) );
        y += yield[i1]*w/norm;
        sumw += w;
      }
    }

    // add second histogram
    for( int i2 = 0; i2 < ptbins2.size()-1 && ptbins2[i2+1]<=1; ++i2 )
    {
      if( ptbins3[i] >= ptbins2[i2] && ptbins3[i+1] <= ptbins2[i2+1] )
      {
        double w = 1.0/square(stat2[i2]/norm2);
        y += yield2[i2]*w/norm2;
        sumw += w;
      }
    }

    y /= sumw;
    yield3.push_back( y );
    error3.push_back( 1.0/sumw );

    const auto binWidth = ptbins3[i+1]-ptbins3[i];
    yield_integral3.push_back( y*binWidth );

    h3->SetBinContent( i+1, yield3[i] );
    h3->SetBinError( i+1, error3[i] );

  }

  // print
  Stream::PrintVector<float>( "float", "pt", &ptbins3[0], ptbins3.size(), "%.3g" );
  Stream::PrintVector<float>( "float", "yield", yield3, "%.3g" );
  Stream::PrintVector<float>( "float", "error", error3, "%.3g" );
  Stream::PrintVector<float>( "float", "yield_int3", yield_integral3, "%.3g" );


  // draw
  {
    auto cv = new TCanvas( "cv", "cv", 800, 800 );
    auto h0 = new TH1F( "dummy", "", 100, 0, 10 );
    h0->GetXaxis()->SetTitle( "p_{#it{T}} (GeV/#it{c})" );
    h0->GetYaxis()->SetTitle( "#it{d}^{2}#it{N}/#it{d}p_{#it{T}}#it{d}#eta" );
    h0->SetMinimum( 3e-8 );
    h0->SetMaximum( 30 );
    h0->Draw();

    h->SetMarkerStyle( 20 );
    h->Draw( "E same" );

    h2->SetMarkerStyle( 20 );
    h2->SetMarkerColor( 2 );
    h2->SetLineColor( 2 );
    h2->Draw( "E same" );

    h3->SetMarkerStyle( 20 );
    h3->SetMarkerColor( 4 );
    h3->SetLineColor( 4 );
    h3->Draw( "E same" );

    gPad->SetLogy( true );
    pdfDocument.Add( cv );
  }
}
