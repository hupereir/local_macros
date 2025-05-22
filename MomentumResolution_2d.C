#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
TString MomentumResolution_2d()
{
  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  const TString tag = "_single_electron";
  const TString inputFile = Form( "DST/CONDOR%s/DST_RECO/dst_reco%s*.root", tag.Data(), tag.Data() );

  static constexpr bool use_micromegas = false;
  const TString pdfFile = use_micromegas ?
    Form( "Figures/MomentumResolution_2d%s_mm.pdf", tag.Data() ):
    Form( "Figures/MomentumResolution_2d%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  // output rootfile
  const TString rootfilename = use_micromegas ?
    Form( "Rootfiles/MomentumResolution_2d%s_mm.root", tag.Data()):
    Form( "Rootfiles/MomentumResolution_2d%s.root", tag.Data());
  RootFile rootFile( rootfilename );

  // variable names
  const TString var( "_tracks._pt/_tracks._truth_pt" );
  const TString var2d = Form( "%s:_tracks._truth_pt", var.Data() );
  TCut cut( "_tracks._nclusters_mvtx>=2 && _tracks._nclusters_intt>=1 && _tracks._nclusters_tpc>=20");
  // TCut cut( "_tracks._nclusters_mvtx>=2 && _tracks._nclusters_intt>=1");
  // TCut cut;
  if( use_micromegas ) cut = cut && "_tracks._nclusters_micromegas>=2";

  static constexpr bool use_fit = true;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // get 2D histogram
//   constexpr int nptbins = 4;
//   const std::array<double, nptbins+1> ptbins = {{ 0.5, 5, 10, 15, 20 }};

  constexpr int nptbins = 19;
  const std::array<double, 20> ptbins = {{ 0.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 }};

//   constexpr int nptbins = 38;
//   const std::array<double, 39> ptbins = {{
//     0.2, 0.3, 0.4, 0.6, 1., 2, 3, 4, 5, 6,
//     7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
//     17, 18, 19, 20, 22, 24, 26, 28, 30, 32,
//     34, 36, 38, 40, 42, 44, 46, 48, 50
//   }};

  // auto h2 = new TH2F( "h2d", "h2d", nptbins, &ptbins[0], 100, 0.8, 1.2 );
  auto h2 = new TH2F( "h2d", "h2d", nptbins, &ptbins[0], 100, 0.7, 1.3 );
  // auto h2 = new TH2F( "h2d", "h2d", nptbins, &ptbins[0], 100, 0.1, 4 );

  // Utils::max_entries = 1e4;

  Utils::TreeToHisto( tree, "h2d", var2d, cut, false );
  h2->SetTitle( "" );
  h2->GetYaxis()->SetTitle( "#it{p}_{T,track}/#it{p}_{T,truth}" );
  h2->GetXaxis()->SetTitle( "#it{p}_{T,truth} (GeV/#it{c})" );
  rootFile.Add( h2 );

  // create TGraph to store resolution vs momentum
  auto tg = new TGraphErrors;
  tg->SetName( "resolution" );

  // create TGraph to store momentum scale vs momentum
  auto tg_mean = new TGraphErrors;
  tg_mean->SetName( "momentum scale" );

  // create canvas
  const TString cvName( "cv" );
  auto cv = new TCanvas( cvName, cvName, 800, 800 );
  Draw::DivideCanvas( cv, nptbins, false );

  // loop over bins
  for( int i = 0; i < nptbins; ++i )
  {

    // get pt slice and fit
    const TString hname = Form( "h_%i", i );
    auto h = h2->ProjectionY( hname, i+1, i+1 );
    std::cout << "MomentumResolution_2d - bin: " << i << " entries: " << h->GetEntries() << std::endl;

    if( !h->GetEntries() ) continue;

    h->Fit( "gaus", "Q0" );
    rootFile.Add( h );

    cv->cd( i+1 );
    h->Draw();
    auto f = h->GetFunction( "gaus" );
    if( !f ) continue;

    f->SetLineColor(2);
    f->Draw( "same" );

    gPad->SetLogy( true );
    Draw::VerticalLine( gPad, 1 )->Draw();

    // add pt range
    Draw::PutText( 0.2, 0.8, Form( "%.1f < #it{p}_{T,truth} < %.1f GeV/#it{c} ", ptbins[i], ptbins[i+1] ) );

    if( use_fit )
    {
      // save momentum scale
      tg_mean->SetPoint( i, 0.5*(ptbins[i]+ptbins[i+1]), f->GetParameter(1));
      tg_mean->SetPointError(i, 0, f->GetParError(1));

      // save resolution in histogram
      tg->SetPoint( i, 0.5*(ptbins[i]+ptbins[i+1]), f->GetParameter(2));
      tg->SetPointError(i, 0, f->GetParError(2));
    } else {
      // save momentum scale
      tg_mean->SetPoint( i, 0.5*(ptbins[i]+ptbins[i+1]), h->GetMean());
      tg_mean->SetPointError(i, 0, h->GetMeanError());

      // save resolution in histogram
      tg->SetPoint( i, 0.5*(ptbins[i]+ptbins[i+1]), h->GetRMS());
      tg->SetPointError(i, 0, h->GetRMSError());
    }
  }
  pdfDocument.Add( cv );

  {
    auto cv = new TCanvas( "cvtgl", "cvtgl", 800, 600 );
    cv->SetLeftMargin( 0.16 );

    auto h = new TH1F( "dummy", "", 100, 0, ptbins[nptbins] );
    h->SetMinimum(0);
    h->SetMaximum(1.2*Utils::GetMaximum( tg ));
    // h->SetMaximum(0.08);
    // h->SetMaximum(0.03);
    // h->SetMaximum(1.2);
    h->GetXaxis()->SetTitle( "#it{p}_{T,truth} (GeV/#it{c})" );
    h->GetYaxis()->SetTitle( "#sigma( #it{p}_{T,track}/#it{p}_{T,truth} )" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();

    tg->SetMarkerStyle(20);
    tg->SetLineColor(1);
    tg->SetMarkerColor(1);
    tg->Draw("P");

    rootFile.Add( tg );
    pdfDocument.Add( cv );

  }

  {
    auto cv = new TCanvas( "cvtgl2", "cvtgl2", 800, 600 );
    cv->SetLeftMargin( 0.16 );

    auto h = new TH1F( "dummy", "", 100, 0, ptbins[nptbins] );
    h->SetMinimum(0.8);
    h->SetMaximum(1.4);
    h->GetXaxis()->SetTitle( "#it{p}_{T,truth} (GeV/#it{c})" );
    h->GetYaxis()->SetTitle( "#LT#it{p}_{T,track}/#it{p}_{T,truth}#GT" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();

    tg_mean->SetMarkerStyle(20);
    tg_mean->SetLineColor(1);
    tg_mean->SetMarkerColor(1);
    tg_mean->Draw("P");


    Draw::HorizontalLine( cv, 1 )->Draw();

    rootFile.Add( tg_mean );
    pdfDocument.Add( cv );

  }

  return pdfFile;
}
