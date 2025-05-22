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
#include <TStyle.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
float delta_phi( float phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

namespace 
{
 
  // pdf document
  PdfDocument pdfDocument;
  
  // initial guess for max residuals
  std::array<float, nDetectors> max_det_residual = {0.003, 0.01, 0.2, 0.2, 0.2, 0.5, 0.5};
  float max_residual = 0;
  
  // variable
  const TString var( "_tracks._clusters[]._trk_r*delta_phi(_tracks._clusters[]._trk_phi - _tracks._clusters[]._truth_phi)" );

  TCut momentum_cut = ("_tracks._pt>0.5" );
  const TCut pattern_cut( 
    "_tracks._truth_pt>0.5"
    "&&_tracks._nclusters_mvtx>0"
    "&&_tracks._nclusters_intt>=2"
    "&&_tracks._nclusters_micromegas==0");

  // save histogram
  std::vector<std::unique_ptr<TH1>> h_array;
    
  //____________________________________________________________________________
  TGraphErrors* get_residuals( TTree* tree, TString tag )
  {
    
    // create TGraph to store resolution vs layer
    auto tg = new TGraphErrors();
    tg->SetName(Form( "residuals_%s", tag.Data()));
    
    // variable names
    const TString var2d = Form( "%s:_tracks._clusters[]._layer", var.Data() );

    // configuration
    const bool do_fit = true;
    
    // loop over detectors
    for( int idet = 0; idet < nDetectors; ++idet )
    {
      
      // skip detector 6, which is z segmented micromegas
      if( idet == 6 ) continue;
      
      const TCut detector_cut( Form( "_tracks._clusters[]._layer>=%i &&_tracks._clusters[]._layer<%i ", firstLayer[idet], firstLayer[idet] + nLayers[idet] ) );
      
      const TString hname( Form( "deltarphi_%i_%s", idet, tag.Data() ) );
      std::unique_ptr<TH2> h2d( new TH2F( hname, "", nLayers[idet], firstLayer[idet], firstLayer[idet] + nLayers[idet], 100, -max_det_residual[idet], max_det_residual[idet] ) );
      Utils::TreeToHisto( tree, hname, var2d, momentum_cut&&pattern_cut&&detector_cut, false );
      
      // create canvas
      const TString cvName = Form( "cv_%i", idet );
      std::unique_ptr<TCanvas> cv( new TCanvas( cvName, cvName, 800, 800 ) );
      Draw::DivideCanvas( cv.get(), nLayers[idet], false );
      
      // loop over layers
      for( int ilayer = 0; ilayer < nLayers[idet]; ++ilayer )
      {
        
        int layerIndex = firstLayer[idet] + ilayer;
        const auto hname = Form( "h_%i_%s", layerIndex, tag.Data() );
        std::unique_ptr<TH1> h( h2d->ProjectionY( hname, ilayer+1, ilayer+1 ) );
        h->SetTitle( hname );
        h->SetLineColor( 1 );
        h->SetMarkerColor( 1 );
        h->GetXaxis()->SetTitle( "r.#Delta#phi_{track-truth} (cm)" );
        h->GetXaxis()->SetRangeUser( -max_det_residual[idet], max_det_residual[idet] );
        
        cv->cd( ilayer+1 );
        h->SetMaximum( 1.2*h->GetMaximum() );
        h->Draw();
        
        // fit
        const auto entries( h->GetEntries() );
        std::cout << "DeltaRPhi_truth - layer: " << layerIndex << " entries: " << entries << std::endl;
        if( entries )
        {
          if( do_fit )
          {
            const auto result = std::min( Fit( h.get() ), Fit_box( h.get() ) );
            if( result._valid )
            {
              auto f = result._function;
              f->Draw("same");
              auto h = f->GetHistogram();
              
              auto mean = f->GetParameter(1);
              auto meanError = f->GetParError(1);
              Draw::PutText( 0.2, 0.8, Form( "mean = %.3g #pm %.3g #mum", mean*1e4, meanError*1e4 ) );
              
              auto rms = h->GetRMS();
              auto error = f->GetParError(2);
              Draw::PutText( 0.2, 0.7, Form( "#sigma = %.3g #pm %.3g #mum", rms*1e4, error*1e4 ) );
              
              tg->SetPoint( layerIndex, radius[layerIndex], rms*1e4 );
              tg->SetPointError( layerIndex, 0, error*1e4 );
              
            } else {
              
              std::cout << "DeltaRPhi_truth - skipping layer " << layerIndex << " (failed fit)" << std::endl;
              
            }
            
          } else {
            
            auto mean = h->GetMean();
            auto meanError = h->GetMeanError();
            Draw::PutText( 0.2, 0.8, Form( "mean = %.3g #pm %.3g #mum", mean*1e4, meanError*1e4 ) );
            
            auto rms = h->GetRMS();
            auto error = h->GetRMSError();
            Draw::PutText( 0.2, 0.7, Form( "#sigma = %.3g #pm %.3g #mum", rms*1e4, error*1e4 ) );
            
            tg->SetPoint( layerIndex, radius[layerIndex], rms*1e4 );
            tg->SetPointError( layerIndex, 0, error*1e4 );
            
          }
        }
        
        // draw vertical line at zero
        gPad->Update();
        Draw::VerticalLine( gPad, 0 )->Draw();
        
        // save in array
        h_array.push_back( std::move(h) );
        
      }
      
      cv->Update();
      cv->cd(0);
      pdfDocument.Add( cv.get() );
      
    }

    return tg;
  
  }
  
}

//____________________________________________________________________________
TString DeltaRPhi_truth_pt( TString tag = TString() )
{

  set_style( false );

//   // if( tag.IsNull() ) tag = "_1k_realistic_full_micromegas_notpc-new" ;
//   const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  if( tag.IsNull() ) tag =  "_flat_full_micromegas_notpc" ;
  
  // need to merge flat and realistic files to have good pT coverage in all 1/pt bins
  std::vector<TString> inputFiles = {
    "DST/CONDOR_realistic_full_micromegas_notpc-new/dst_eval_realistic_full_micromegas_notpc_*.root",
    "DST/CONDOR_flat_full_micromegas_notpc-new/dst_eval_flat_full_micromegas_notpc_*.root"
  };

  const TString pdfFile = Form( "Figures/DeltaRPhi_truth%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/DeltaRPhi_truth%s.root", tag.Data() );

  std::cout << "DeltaRPhi_truth - pdfFile: " << pdfFile << std::endl;

  pdfDocument = PdfDocument( pdfFile );

  // file manager
  FileManager fileManager;
  for( const auto& inputFile:inputFiles )
  { 
    std::cout << "DeltaRPhi_truth - inputFile: " << inputFile << std::endl;
    fileManager.AddFiles( inputFile ); 
  }
  
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  std::cout << "DeltaRPhi_truth_pt - finished loading tree" << std::endl;
  
  // optimize max residual
  for( int idet = 0; idet < nDetectors; ++idet )
  {
    
    // skip detector 6, which is z segmented micromegas
    if( idet == 6 ) continue;

    const TString hname( Form( "deltarphi_%i_0", idet ) );
    const TCut layer_cut( Form( "_tracks._clusters[]._layer==%i", firstLayer[idet] ) );

    for( int i=0; i<3; ++i )
    {
      std::unique_ptr<TH1> h1( new TH1F( hname, "", 500, -max_det_residual[idet], max_det_residual[idet] ) );
      Utils::TreeToHisto( tree, hname, var, momentum_cut&&layer_cut&&pattern_cut, false );
      max_det_residual[idet] = 8*h1->GetRMS() + std::abs( h1->GetMean() );;
    }

  }
  
  // max residual
  max_residual = *std::max_element( max_det_residual.cbegin(), max_det_residual.cend() )/5;
  std::cout << "DeltaRPhi_truth_pt - finished setting residuals" << std::endl;

  std::unique_ptr<TCanvas> cv( new TCanvas( "cvtg", "cvtg", 800, 600 ) );
  cv->SetLeftMargin( 0.16 );
  
  // std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 0, 90 ) );
  std::unique_ptr<TH1> h( new TH1F( "dummy", "", 100, 20, 80 ) );
  h->SetMinimum(0);
  h->SetMaximum( 600 );
  h->GetXaxis()->SetTitle( "r (cm)" );
  h->GetYaxis()->SetTitle( "#sigma_{r.#Delta#phi} (track-truth) (#mum)" );
  h->GetYaxis()->SetTitleOffset( 1.6 );
  h->Draw();

  constexpr int nptbins = 20;
  std::array<float, nptbins+1> invptbins = {{ 
    0, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 
    1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05 }};
  constexpr std::array<int, nptbins> color = {{ 
    kBlue, kCyan+2, kGreen+2, kOrange+1, kRed, kBlue, kCyan+2, kGreen+2, kOrange+1, kRed, 
    kBlue, kCyan+2, kGreen+2, kOrange+1, kRed, kBlue, kCyan+2, kGreen+2, kOrange+1, kRed
  }};

  // constexpr int nptbins = 2;
  // std::array<float, nptbins+1> invptbins = {{ 0, 0.5, 1.0 }};
  // constexpr std::array<int, nptbins> color = {{ kBlue, kCyan+2 }};

  // constexpr int nptbins = 1;
  // std::array<float, nptbins+1> invptbins = {{ 0, 2.0 }};
  // constexpr std::array<int, nptbins> color = {{ kBlue }};
    
  std::vector<TGraphErrors*> tg_array;  
  for( int ibin = 0; ibin < nptbins; ++ibin )
  {
    std::cout << "DeltaRPhi_truth_pt.C - range: " << invptbins[ibin] << ", " << invptbins[ibin+1] << std::endl;
    momentum_cut = Form( "_tracks._pt>0.5 && (1./_tracks._truth_pt) >= %f && (1./_tracks._truth_pt) < %f", 
      invptbins[ibin], 
      invptbins[ibin+1] );
  
    auto tg = get_residuals( tree, Form( "%i", ibin ));
    
    cv->cd();
    tg->SetMarkerStyle(20);
    tg->SetLineColor(color[ibin]);
    tg->SetMarkerColor(color[ibin]);
    tg->Draw("P");
    
    tg_array.push_back( tg );
        
  }

  // store to pdf document
  pdfDocument.Add( cv.get() );

  // save everything in rootfiles
  std::unique_ptr<TFile> output( TFile::Open( rootFile, "RECREATE" ) );
  output->cd();
  for( auto&& h:h_array) { if(h) h->Write(); }
  for( auto&& tg:tg_array) { if(tg) tg->Write(); }
  output->Close();
  
  h_array.clear();
  
  return pdfFile;

}
