#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TChain.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
TString DeltaZ_cluster_tpc()
{

  set_style( false );

  // initial guess for max residuals
  // std::array<float, nDetectors> max_det_residual = { 0.003, 5, 1.2, 1.2, 1.2, 0.1, 0.1};
  std::array<float, nDetectors> max_det_residual = { 0.003, 5, 0.5, 0.5, 0.5, 0.1, 0.1};

  // input files
  const TString tag = "_genfit_truth_notpc_nodistortion-new";
  const TString inputFile = Form( "DST/CONDOR_realistic_micromegas/dst_reco%s/dst_reco*_?.root", tag.Data() );
  
  const TString pdfFile = Form( "Figures/DeltaZ_cluster_tpc%s.pdf", tag.Data() );

  std::cout << "DeltaZ_cluster - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaZ_cluster - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // configuration
  const bool do_fit = true;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  // variable names
  const TString var( "_tracks._clusters._z - _tracks._clusters._truth_z" );
  const TString var2d = Form( "%s:_tracks._clusters._layer", var.Data() );

  const TCut momentum_cut;
  
  static constexpr int ncuts = 2;
  const std::array<TCut, ncuts> cluster_cut = {{
    TCut( "_tracks._clusters._z > 0" ),
    TCut( "_tracks._clusters._z < 0" )
  }};
  
  const std::array<TString,ncuts> cut_tag = {{ "negz", "posz" }};
  const std::array<TString,ncuts> label = {{
    "z_{clus} > 0",
    "z_{clus} < 0",
  }};
  
  const std::array<int, ncuts> color = {{ 1, 2 }};
    
  // TGraphs
  std::array<TGraphErrors*, ncuts> tg;
  std::array<TGraphErrors*, ncuts> tgl;
  std::array<TGraphErrors*, ncuts> tgm;
  std::array<TGraphErrors*, ncuts> tgml;

  for( int icut = 0; icut < ncuts; ++icut )
  {
    tg[icut] = new TGraphErrors();
    tg[icut]->SetName( Form( "residuals %s", cut_tag[icut].Data() ) );
    tg[icut]->SetMarkerStyle(20);
    tg[icut]->SetLineColor(color[icut]);
    tg[icut]->SetMarkerColor(color[icut]);

    tgl[icut] = new TGraphErrors();
    tgl[icut]->SetName( Form( "residuals_layers %s", cut_tag[icut].Data() ) );
    tgl[icut]->SetMarkerStyle(20);
    tgl[icut]->SetLineColor(color[icut]);
    tgl[icut]->SetMarkerColor(color[icut]);

    tgm[icut] = new TGraphErrors();
    tgm[icut]->SetName( Form( "residuals_mean %s", cut_tag[icut].Data() ) );
    tgm[icut]->SetMarkerStyle(20);
    tgm[icut]->SetLineColor(color[icut]);
    tgm[icut]->SetMarkerColor(color[icut]);

    tgml[icut] = new TGraphErrors();
    tgml[icut]->SetName( Form( "residuals_mean_layer %s", cut_tag[icut].Data() ) );
    tgml[icut]->SetMarkerStyle(20);
    tgml[icut]->SetLineColor(color[icut]);
    tgml[icut]->SetMarkerColor(color[icut]);
  };
  
//   // optimize max residual
//   for( int idet = 0; idet < nDetectors; ++idet )
//   {
// 
//     // skip detector 5, which is phi segmented micromegas
//     if( idet == 5 ) continue;
// 
//     const TString hname( Form( "deltarphi_%i_0", idet ) );
//     const TCut layer_cut( Form( "_tracks._clusters._layer==%i", firstLayer[idet] ) );
// 
//     for( int i=0; i<3; ++i )
//     {
//       std::unique_ptr<TH1> h1( new TH1F( hname, "", 500, -max_det_residual[idet], max_det_residual[idet] ) );
//       Utils::TreeToHisto( tree, hname, var, momentum_cut&&layer_cut, false );
//       max_det_residual[idet] = 5*h1->GetRMS() + std::abs(h1->GetMean() );
//     }
// 
//   }

  // loop over detectors
  for( int idet = 0; idet < nDetectors; ++idet )
  {
    
    std::cout << "DeltaZ_cluster - detector " << idet << " max_det_residual: " << max_det_residual[idet] << std::endl;
    
    // skip detector 5, which is phi segmented micromegas
    if( idet == 5 ) continue;

    const TCut layer_cut = Form( "_tracks._clusters._layer >= %i && _tracks._clusters._layer < %i", firstLayer[idet], firstLayer[idet] + nLayers[idet] );

    std::array<std::unique_ptr<TH2>,ncuts> h2d;
    for( int icut = 0; icut < ncuts; ++icut )
    {
      const TString hname( Form( "deltaz_%i_%s", idet, cut_tag[icut].Data() ) );
      h2d[icut].reset( new TH2F( hname, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 100, -max_det_residual[idet], max_det_residual[idet] ) );
      Utils::TreeToHisto( tree, hname, var2d, cluster_cut[icut]&&momentum_cut&&layer_cut, false );
    }
    
    // create canvas
    const TString cvName = Form( "cv_%i", idet );
    std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, (idet == 0 ) ? 400:800 ) );
    Draw::DivideCanvas( cv.get(), nLayers[idet], false );

    // loop over layers
    for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
    {

      int layerIndex = firstLayer[idet] + ilayer;
      
      bool first = true;
      for( int icut = 0; icut < ncuts; ++icut )
      {
        const auto hname = Form( "h_%i_%i", layerIndex,icut );
        // std::unique_ptr<TH1> h( h2d[icut]->ProjectionY( hname, ilayer+1, ilayer+1 ) );
        auto h( h2d[icut]->ProjectionY( hname, ilayer+1, ilayer+1 ) );
        if( !h->GetEntries() ) continue;
        
        h->SetTitle( hname );
        h->SetLineColor( color[icut] );
        h->SetMarkerColor( color[icut] );
        h->GetXaxis()->SetTitle( "#Deltaz_{cluster-truth} (cm)" );
        h->GetXaxis()->SetRangeUser( -max_det_residual[idet], max_det_residual[idet] );
        h->SetMaximum( h->GetMaximum()*1.2 );

        cv->cd( ilayer+1 );
        if( first ) { first = false; h->Draw(); }
        else h->Draw( "same" );

        if( idet == 0 && icut == 0 )
        {
          h->GetXaxis()->SetMaxDigits( 2 );
          gPad->SetRightMargin( 0.12 );
        }

        // fit
        const auto entries( h->GetEntries() );
        std::cout << "DeltaRPhi_cluster - layer: " << layerIndex << " entries: " << entries << std::endl;
        if( entries )
        {
          if( do_fit )
          {
            
            const auto result = std::min( Fit(h), Fit_box(h) );
            
            if( result._valid )
            {
              auto f = result._function;
              f->SetLineColor( color[icut] );
              f->Draw("same");
              auto h = f->GetHistogram();
              
              auto mean = f->GetParameter(1);
              auto meanError = f->GetParError(1);
              
              tgml[icut]->SetPoint( layerIndex, layerIndex, mean );
              tgml[icut]->SetPointError( layerIndex, 0, meanError );
              
              tgm[icut]->SetPoint( layerIndex, radius[layerIndex], mean );
              tgm[icut]->SetPointError( layerIndex, 0, meanError );
              
              auto rms = h->GetRMS();
              auto error = f->GetParError(2);
              Draw::PutText( 0.2, 0.8 - 0.1*icut, Form( "%s, mean = %.3g #pm %.3g #mum, #sigma = %.3g #pm %.3g #mum", label[icut].Data(), mean*1e4, meanError*1e4, rms*1e4, error*1e4 ) );

              // add to histogram
              
              tgl[icut]->SetPoint( layerIndex, layerIndex, rms*1e4 );
              tgl[icut]->SetPointError( layerIndex, 0, error*1e4 );

              tg[icut]->SetPoint( layerIndex, radius[layerIndex], rms*1e4 );
              tg[icut]->SetPointError( layerIndex, 0, error*1e4 );
              
            } else {
            
              std::cout << "DeltaZ - skipping layer " << layerIndex << " (failed fit)" << std::endl;
            
            }
            
          } else {
            
            auto mean = h->GetMean();
            auto meanError = h->GetMeanError();
            
            tgml[icut]->SetPoint( layerIndex, layerIndex, mean );
            tgml[icut]->SetPointError( layerIndex, 0, meanError );
            
            tgm[icut]->SetPoint( layerIndex, radius[layerIndex], mean );
            tgm[icut]->SetPointError( layerIndex, 0, meanError );
            
            auto rms = h->GetRMS();
            auto error = h->GetRMSError();
            Draw::PutText( 0.2, 0.8 - 0.1*icut, Form( "%s, mean = %.3g #pm %.3g #mum, #sigma = %.3g #pm %.3g #mum", label[icut].Data(), mean*1e4, meanError*1e4, rms*1e4, error*1e4 ) );
            
            tgl[icut]->SetPoint( layerIndex, layerIndex, rms*1e4 );
            tgl[icut]->SetPointError( layerIndex, 0, error*1e4 );
            
            tg[icut]->SetPoint( layerIndex, radius[layerIndex], rms*1e4 );
            tg[icut]->SetPointError( layerIndex, 0, error*1e4 );
          }
        }
        
      }

      // draw vertical line at zero
      gPad->Update();
      Draw::VerticalLine( gPad, 0 )->Draw();
    }

    cv->Update();
    cv->cd(0);
    pdfDocument.Add( cv.get() );

  }

  // RMS tgraph
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtgl", "cvtgl", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, nLayersTotal ) );
    h->SetMinimum(0);
    h->SetMaximum(1.2*Utils::GetMaximum(tgl[0]));
    h->GetXaxis()->SetTitle( "layer id" );
    h->GetYaxis()->SetTitle( "#sigma_{z} (cluster-truth) (#mum)" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();

    for( int icut = 0; icut < ncuts; ++icut )
    { tgl[icut]->Draw("P"); }

    pdfDocument.Add( cv.get() );
  }

  // RMS tgraph
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtg", "cvtg", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, 90 ) );
    h->SetMinimum(0);
    h->SetMaximum(1.2*Utils::GetMaximum(tg[0]));
    h->GetXaxis()->SetTitle( "r (cm)" );
    h->GetYaxis()->SetTitle( "#sigma_{z} (cluster-truth) (#mum)" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();

    for( int icut = 0; icut < ncuts; ++icut )
    { tg[icut]->Draw("P"); }

    pdfDocument.Add( cv.get() );
  }

  // static constexpr float max_offset = 1.;
  static constexpr float max_offset = 0.04;

  // mean tgraph
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtgl", "cvtgl", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, nLayersTotal ) );
    h->SetMinimum(-max_offset);
    h->SetMaximum(max_offset);
    h->GetXaxis()->SetTitle( "layer id" );
    h->GetYaxis()->SetTitle( "#LT#Deltaz#GT (cluster-truth) (#mum)" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();
    
    auto legend = new TLegend( 0.2, 0.82, 0.60, 0.93, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

    for( int icut = 0; icut < ncuts; ++icut )
    { 
      tgml[icut]->Draw("P"); 
      tgml[icut]->Fit( "pol0", "", "", 7, 54 );
      legend->AddEntry( tgml[icut], label[icut], "PL" );
      
    }

    
    
    pdfDocument.Add( cv.get() );
  }

  // Mean tgraph
  {
    std::unique_ptr<TCanvas> cv( new TCanvas( "cvtg", "cvtg", 800, 600 ) );
    cv->SetLeftMargin( 0.16 );

    std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, 90 ) );
    h->SetMinimum(-max_offset);
    h->SetMaximum(max_offset);
    h->GetXaxis()->SetTitle( "r (cm)" );
    h->GetYaxis()->SetTitle( "#LT#Deltaz#GT (cluster-truth) (#mum)" );
    h->GetYaxis()->SetTitleOffset( 1.6 );
    h->Draw();
   
    auto legend = new TLegend( 0.2, 0.82, 0.60, 0.93, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

    for( int icut = 0; icut < ncuts; ++icut )
    { 
      tgm[icut]->Draw("P");   
      legend->AddEntry( tgm[icut], label[icut], "PL" );
    }

    pdfDocument.Add( cv.get() );
  }
  
  return pdfFile;

}
