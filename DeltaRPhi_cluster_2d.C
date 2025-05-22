#include <RootUtil/FileManager.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/PdfDocument.h>
#include <g4eval/TrackingEvaluator_hp.h>

#include <TTree.h>

R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

namespace
{
  
  template< class T>
    inline constexpr T square( const T& x ) { return x*x; }
    
  template< class T>
    inline T delta_phi( const T& phi )
  {
    if( phi >= M_PI ) return phi - 2*M_PI;
    else if( phi < -M_PI ) return phi + 2*M_PI;
    else return phi;
  }
  
}

TString DeltaRPhi_cluster_2d()
{
 
  // maximum residual
  // const double max_residual = 0.5;
  const double max_residual = 4;
  
  auto hentries = new TH2F( "entries", "entries", 100, -105.5, 105.5, nLayers_tpc, &tpc_radius[0] );
  auto hmean = new TH2F( "mean", "mean", 100, -105.5, 105.5, nLayers_tpc, &tpc_radius[0] );
  auto hrms = new TH2F( "rms", "rms", 100, -105.5, 105.5, nLayers_tpc, &tpc_radius[0] );

  for( const auto& h:{hentries, hmean, hrms} )
  { 
    h->GetXaxis()->SetTitle( "z (cm)" );
    h->GetYaxis()->SetTitle( "r (cm)" );
  }
  
  hentries->GetZaxis()->SetTitle( "entries" );
  hmean->GetZaxis()->SetTitle( "#LTr#Delta#phi#GT" );
  hrms->GetZaxis()->SetTitle( "#sigma_{#Delta#phi}" );
      
//   // input files
//   const TString tag = "_realistic_truth_genfit";
//   const TString inputFile = Form( "DST/dst_reco%s.root", tag.Data() );

  // input files
  const TString tag = "_directlasers-nominal";
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  TString rootfilename = Form( "Rootfiles/DeltaRPhi_cluster_2d%s.root", tag.Data() );
  RootFile rootFile( rootfilename );

  TString pdffilename = Form( "Figures/DeltaRPhi_cluster_2d%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdffilename );
  
  // get tree
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return pdffilename;

  auto container = new TrackingEvaluator_hp::Container;
  tree->SetBranchAddress( "DST#EVAL#TrackingEvaluator_hp::Container", &container );
  
  int total_clusters = 0;

  // loop over entries
  // const int entries = 5e5;
  const int entries = tree->GetEntries();
  std::cout << "DeltaRPhi_cluster_2d - entries: " << entries << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    tree->GetEntry(i);
    if( !(i%100) )
    { std::cout << "DeltaRPhi_cluster_2d - entry: " << i << std::endl; }

    // loop over clusters
    for( const auto& cluster:container->clusters() )
    {
      
      // keep only tpc clusters
      if( cluster._layer < 7 || cluster._layer >= 55 ) continue;
      
      // cluster position
      const auto r = cluster._r;
      const auto z = cluster._z;
      
      // get delta rphi
      const auto drp = cluster._r*( delta_phi( cluster._phi - cluster._truth_phi ) );
      if( std::abs( drp ) > max_residual ) continue;
      
      // fill histograms
      hentries->Fill( z, r );
      hmean->Fill( z, r, drp );
      hrms->Fill( z, r, square(drp) );
      
      ++total_clusters;
    }
    
  }
  
  // finish histogram normalization
  for( int iz = 0; iz < hentries->GetNbinsX(); ++iz )
    for( int ir = 0; ir < hentries->GetNbinsY(); ++ir )
  {
    const auto entries = hentries->GetBinContent( iz+1, ir+1 );
    if( entries > 10 )
    {
      
      const auto mean = hmean->GetBinContent( iz+1, ir+1 )/entries;
      const auto rms = std::sqrt( hrms->GetBinContent( iz+1, ir+1 )/entries - square( mean ) );
      
      // reassign
      hmean->SetBinContent( iz+1, ir+1, mean );
      hrms->SetBinContent( iz+1, ir+1, rms );
    
    } else { 
      
      hentries->SetBinContent( iz+1, ir+1, 0 );
      hmean->SetBinContent( iz+1, ir+1, 0 );
      hrms->SetBinContent( iz+1, ir+1, 0 );
      
    }
    
  }
  
  // save to root files
  rootFile.Add( hentries );
  rootFile.Add( hmean );
  rootFile.Add( hrms );
  
  int i = 0;
  for( const auto& h:{hentries, hmean, hrms })
  {
    const auto cvname = Form( "cv%i", i++ );
    auto cv = new TCanvas( cvname, cvname, 800, 800 );
    cv->SetRightMargin(0.2);
    h->Draw("colz" );
    pdfDocument.Add(cv);
  }

  std::cout << "DeltaRPhi_cluster_2d - cluster statistics total: " << total_clusters << std::endl;

  return pdffilename;

}
