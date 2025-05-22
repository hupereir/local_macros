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
int get_charge( int pid )
{
  
  switch( pid ) 
  { 
    case 211: return 1;
    case -211: return -1;
    default: return 0;
  }
}

//____________________________________________________________________________
TString MomentumResolution_2d_charge( TString tag = TString() )
{
  set_style( false );
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // if( tag.IsNull() ) tag = "_flat_truth_micromegas_corrected_all-coarse";
  if( tag.IsNull() ) tag = "_flat_truth_micromegas_corrected_mm-coarse_extrapolated";
  const TString inputFile = Form( "DST/CONDOR%s/dst_eval*.root", tag.Data() );

  static constexpr bool use_micromegas = false;
  const TString pdfFile = use_micromegas ?
    Form( "Figures/MomentumResolution_2d_charge%s_mm.pdf", tag.Data() ):
    Form( "Figures/MomentumResolution_2d_charge%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  RootFile rootFile( use_micromegas ?
    Form( "Rootfiles/MomentumResolution_2d_charge%s_mm.root", tag.Data() ):
    Form( "Rootfiles/MomentumResolution_2d_charge%s.root", tag.Data() ) );

  // variable names
  const TString var( "_tracks._pt/_tracks._truth_pt" );
  const TString var2d = Form( "%s:_tracks._truth_pt", var.Data() );
  const TCut cut = use_micromegas ?
    TCut( "abs(_tracks._truth_pt)>0.5&&_tracks._nclusters_micromegas>=2" ):
    TCut( "abs(_tracks._truth_pt)>0.5" );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // get 2D histogram
//   constexpr int nptbins = 4;
//   const std::array<double, nptbins+1> ptbins = {{ 0.5, 5, 10, 15, 20 }};
  
  constexpr int nptbins = 19;
  const std::array<double, 20> ptbins = {{ 0.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 }};
  // auto h2 = new TH2F( "h2d", "h2d", nptbins, &ptbins[0], 100, 0.9, 1.1 );
  
  const int ncharge_cut = 2;
  const std::array<TCut,2> charge_cut = {{ "_tracks._charge < 0", "_tracks._charge>0" }};
  // const std::array<TCut,2> charge_cut = {{ "get_charge(_tracks._pid) < 0", "get_charge(_tracks._pid)>0" }};
  const std::array<int,2> color = {{1,2}};
  
  std::array<TH2*, 2> h2 = 
  {{
//     new TH2F( "h2dm", "h2dm", nptbins, &ptbins[0], 100, 0., 2.5 ),
//     new TH2F( "h2dp", "h2dp", nptbins, &ptbins[0], 100, 0., 2.5 )
    new TH2F( "h2dm", "h2dm", nptbins, &ptbins[0], 100, 0.5, 2 ),
    new TH2F( "h2dp", "h2dp", nptbins, &ptbins[0], 100, 0.5, 2 )
  }};

  std::array<TString, 2> label = { "q<0", "q>0" };
  
  std::array<TGraphErrors*, 2> tg = 
  {{
    new TGraphErrors(),
    new TGraphErrors()
  }};

  tg[0]->SetName( "ResolutionM" );
  tg[1]->SetName( "ResolutionP" );

  std::array<TGraphErrors*, 2> tg_mean = 
  {{
    new TGraphErrors(),
    new TGraphErrors()
  }};

  tg_mean[0]->SetName( "MeanM" );
  tg_mean[1]->SetName( "MeanP" );  
  
  // project
  for( int i = 0; i < ncharge_cut; ++i )
  {
    Utils::TreeToHisto( tree, h2[i]->GetName(), var2d, cut && charge_cut[i], false );
    h2[i]->SetTitle( "" );
    h2[i]->GetYaxis()->SetTitle( "#it{p}_{T,track}/#it{p}_{T,truth}" );
    h2[i]->GetXaxis()->SetTitle( "#it{p}_{T,truth} (GeV/#it{c})" );
    rootFile.Add( h2[i] );
  }
    
  // create TGraph to store resolution vs momentum

  // create canvas
  const TString cvName( "cv" );
  auto cv = new TCanvas( cvName, cvName, 800, 800 );
  Draw::DivideCanvas( cv, nptbins, false );
  
  auto legend = new TLegend( 0.5, 0.84, 0.97, 0.94, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  
  // loop over bins
  for( int i = 0; i < nptbins; ++i )
  {
    
    // lover over charge cut
    for( int j = 0; j<ncharge_cut; ++j )
    {
    
      // get pt slice and fit
      const auto hname = Form( "h_%i_%i", i, j );
      auto h = h2[j]->ProjectionY( hname, i+1, i+1 );
      std::cout << "MomentumResolution_2d - bin: " << i << " entries: " << h->GetEntries() << std::endl;
      h->Fit( "gaus", "Q0" );
      h->SetMaximum( 1.2*h->GetMaximum() );
      
      // store
      rootFile.Add( h );
    
      cv->cd( i+1 );
      h->SetLineColor( color[j] );
      if( j ) h->Draw( "same" );
      else h->Draw();
      
      auto f = h->GetFunction( "gaus" );
      f->SetLineColor( color[j] );
      f->Draw( "same" );
      
      if( !j )
      {
        gPad->SetLogy( true );
        Draw::VerticalLine( gPad, 1 )->Draw();
      
        // add pt range
        Draw::PutText( 0.2, 0.8, Form( "%.1f < #it{p}_{T,truth} < %.1f GeV/#it{c} ", ptbins[i], ptbins[i+1] ) );
      }
      
      // save resolution in tgraph
      tg[j]->SetPoint( i, 0.5*(ptbins[i]+ptbins[i+1]), f->GetParameter(2));
      tg[j]->SetPointError(i, 0, f->GetParError(2));    

      // save resolution in tgraph
      tg_mean[j]->SetPoint( i, 0.5*(ptbins[i]+ptbins[i+1]), f->GetParameter(1));
      tg_mean[j]->SetPointError(i, 0, f->GetParError(1));    
      
      // update legend
      if( i == 0 )
      { 
        legend->AddEntry( h, label[j], "L" );
        if( j == 1 ) legend->Draw();
      }
    
    }
    
  }
  pdfDocument.Add( cv );

  {
    auto cv = new TCanvas( "cvtgl", "cvtgl", 600, 600 );
    cv->SetLeftMargin( 0.16 );
    
    auto h = new TH1F( "dummy", "", 100, 0, ptbins[nptbins] );
    h->SetMinimum(0);
    h->SetMaximum(1.2*Utils::GetMaximum( tg[0] ));
    h->GetXaxis()->SetTitle( "#it{p}_{T,truth} (GeV/#it{c})" );
    h->GetYaxis()->SetTitle( "#sigma( #it{p}_{T,track}/#it{p}_{T,truth} )" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();
    
    for( int i = 0; i < ncharge_cut; ++i )
    {
      tg[i]->SetMarkerStyle(20);
      tg[i]->SetLineColor(color[i]);
      tg[i]->SetMarkerColor(color[i]);
      tg[i]->Draw("P");

      rootFile.Add( tg[i] );
      
    }
    
    legend->Draw();
    pdfDocument.Add( cv );
    
  }
  
  {
    auto cv = new TCanvas( "cvtg2", "cvtg2", 600, 600 );
    cv->SetLeftMargin( 0.16 );
    
    auto h = new TH1F( "dummy", "", 100, 0, ptbins[nptbins] );
    h->SetMinimum(0);
    h->SetMaximum(2.5);
    h->GetXaxis()->SetTitle( "#it{p}_{T,truth} (GeV/#it{c})" );
    h->GetYaxis()->SetTitle( "< #it{p}_{T,track}/#it{p}_{T,truth} >" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();
    
    for( int i = 0; i < ncharge_cut; ++i )
    {
      tg_mean[i]->SetMarkerStyle(20);
      tg_mean[i]->SetLineColor(color[i]);
      tg_mean[i]->SetMarkerColor(color[i]);
      tg_mean[i]->Draw("P");

      rootFile.Add( tg_mean[i] );
      
    }
    
    Draw::HorizontalLine( cv, 1 )->Draw();
    
    legend->Draw();
    pdfDocument.Add( cv );
    
  }

  return pdfFile;
}
