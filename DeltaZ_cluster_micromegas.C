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

#include "LayerDefines.h"
#include "Fit.C"

R__LOAD_LIBRARY(libRootUtilBase.so)

//____________________________________________________________________________
TString DeltaZ_cluster_micromegas()
{

  set_style( false );

  // input files
  const TString tag = "_acts_truth_nodistortion";
  const TString inputFile = Form( "DST/CONDOR_realistic_micromegas/dst_reco%s/dst_reco*_??.root", tag.Data() );

//   // input files
//   const TString tag = "_truth_no_distortion";
//   const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaZ_cluster_micromegas%s.pdf", tag.Data() );
  const TString rootFile  = Form( "Rootfiles/DeltaZ_cluster%s.root", tag.Data() );

  std::cout << "DeltaZ_cluster - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaZ_cluster - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  if( false )
  {
    
    // variable names
    // #define USE_TRACKS
    #ifdef USE_TRACKS
    const TString var( "_tracks._clusters._z - _tracks._clusters._truth_z" );
    const TString var2d = Form( "%s:_tracks._clusters._z_size", var.Data() );
    #else
    const TString var( "_clusters._z - _clusters._truth_z" );
    const TString var2d = Form( "%s:_clusters._z_size", var.Data() );
    #endif
    
    // const TCut cluster_cut( "_clusters._truth_size == 3" );
    const TCut cluster_cut;
    const TCut momentum_cut;
    
    const int layer = 56;

    // layer cut
    #ifdef USE_TRACKS
    const TCut layer_cut = Form( "_tracks._clusters._layer == %i", layer );
    #else
    const TCut layer_cut = Form( "_clusters._layer == %i", layer );
    #endif
  
    const int csize_min = 1;
    const int csize_max = 6;
    const int n_csize = csize_max - csize_min;
    
    auto cv = new TCanvas( "cv", "cv", 800, 800 );
    Draw::DivideCanvas( cv, n_csize );
    
    const TString hname = "delta_z_2d";
    auto h2d = new TH2F( hname, "", n_csize, csize_min, csize_max, 100, -0.2, 0.2 );
    Utils::TreeToHisto( tree, hname, var2d, cluster_cut&&momentum_cut&&layer_cut, false );
    h2d->GetXaxis()->SetTitle( "cluster size" );
    h2d->GetYaxis()->SetTitle( "#Deltaz_{clus-truth} (cm)" );

    std::cout << "DeltaZ_cluster_micromegas - entries: " << h2d->GetEntries() << std::endl;
        
    // also store RMS and error in separate TGraph
    auto tg = new TGraphErrors;

    // project to cluster size
    for( int csize = csize_min; csize < csize_max; ++csize )
    {
            
      int icsize = csize-csize_min;
      const auto hname = Form( "h_%i", icsize );
      auto h( h2d->ProjectionY( hname, icsize+1, icsize+1 ) );
      
      std::cout << "DeltaZ_cluster_micromegas -"
        << " csize: " << csize 
        << " entries: " << h->GetEntries()
        << std::endl;

      cv->cd( icsize+1 );
      h->Draw();
      
      gPad->Update();
      Draw::PutText( 0.2, 0.85, Form( "cluster size: %i", csize ) );
      Draw::PutText( 0.2, 0.75, Form( "entries: %.0f", h->GetEntries() ) );
      Draw::PutText( 0.2, 0.65, Form( "RMS: %.5f", h->GetRMS() ) );
           
      tg->SetPoint( icsize, csize, h->GetRMS() );
      tg->SetPointError( icsize, 0, h->GetRMSError() );

    }
    
    pdfDocument.Add( cv );

    {
      // draw tgraph
      auto cv = new TCanvas( "cv1", "cv1", 800, 800 );
      cv->SetLeftMargin(0.2);
      
      tg->GetXaxis()->SetTitle( "cluster size" );
      tg->GetYaxis()->SetTitle( "RMS(#Deltaz_{clus-truth}) (cm)" );
      tg->GetYaxis()->SetTitleOffset(2.0);
      tg->SetMinimum(0);
      tg->SetMaximum(0.07);
      tg->SetMarkerStyle(20);
      tg->Draw("AP" );
     
      gPad->Update();
      
      pdfDocument.Add( cv );
    }
    
  }

  if( true )
  {
    
    // residuals vs track angle
    #ifdef USE_TRACKS
    const TString var( "_tracks._clusters._z - _tracks._clusters._truth_z" );
    const TString var2d = Form( "%s:fabs(_tracks._clusters._truth_beta)", var.Data() );
    #else
    const TString var( "_clusters._z - _clusters._truth_z" );
    const TString var2d = Form( "%s:fabs(_clusters._truth_beta)", var.Data() );
    #endif
    
    // const TCut cluster_cut( "_clusters._truth_size == 3" );
    const TCut cluster_cut;
    const TCut momentum_cut;
    const int layer = 56;

    // layer cut
    #ifdef USE_TRACKS
    const TCut layer_cut = Form( "_tracks._clusters._layer == %i", layer );
    #else
    const TCut layer_cut = Form( "_clusters._layer == %i", layer );
    #endif
  
    auto cv = new TCanvas( "cv2", "cv2", 800, 800 );
    cv->SetLeftMargin(0.2);
    
    const TString hname = "delta_rphi_beta_2d";
    const int nanglebins = 20;
    const double minbeta = 0;
    const double maxbeta = 0.85;
    auto h2d = new TH2F( hname, "", nanglebins, minbeta, maxbeta, 100, -0.1, 0.1 );
    Utils::TreeToHisto( tree, hname, var2d, cluster_cut&&momentum_cut&&layer_cut, false );
    h2d->GetXaxis()->SetTitle( "tan(#beta)" );
    h2d->GetYaxis()->SetTitle( "r.#Delta#phi_{clus-truth} (cm)" );
    h2d->GetYaxis()->SetTitleOffset(2.0);

    const TString pname = "delta_rphi_beta_p";
    auto p = new TProfile( pname, "", nanglebins, minbeta, maxbeta );
    p->SetErrorOption( "s" );
    Utils::TreeToHisto( tree, pname, var2d, cluster_cut&&momentum_cut&&layer_cut, false );
    p->SetLineColor( 2 );
    
    h2d->Draw();
    p->Draw("same");
    cv->Update();
    pdfDocument.Add( cv );

    // store RMS vs track angle
    auto tg = new TGraphErrors;
    for( int i = 0; i < p->GetNbinsX(); ++i )
    { tg->SetPoint( i, p->GetXaxis()->GetBinCenter(i+1), p->GetBinError( i+1 ) ); }
 
    {
      // draw tgraph
      auto cv = new TCanvas( "cv1", "cv1", 800, 800 );
      cv->SetLeftMargin(0.2);
      
      tg->GetXaxis()->SetTitle( "tan(#alpha)" );
      tg->GetYaxis()->SetTitle( "RMS(#Deltaz_{clus-truth}) (cm)" );
      tg->GetYaxis()->SetTitleOffset(2.0);
      tg->SetMinimum(0);
      tg->SetMaximum(0.07);
      tg->SetMarkerStyle(20);
      tg->Draw("AP" );
      
//       tg->Fit("pol1", "0Q");
//       auto f = tg->GetFunction("pol1");
//       f->SetLineColor(2);
//       f->Draw("same");
     
      gPad->Update();
      
      pdfDocument.Add( cv );
    }

  }
  
  return pdfFile;
}
