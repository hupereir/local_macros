#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>

namespace
{
 
  // detector names
  std::vector<std::string> detector_names;

  MicromegasMapping mapping;
 
  // save detector names
  void save_detector_names()
  {
    
    // get detector names that match tile/layer
    for( int ilayer = 0; ilayer < 2; ++ilayer )
      for( int itile = 0; itile < 8; ++itile )
    {
      const int layer = 55+ilayer;
      const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
      detector_names.push_back( std::move( name ) );
    }
  }

  int get_color( int layer, int tile, int region )
  {
    std::array<int, 8> color_base = {kOrange, kRed, kPink, kMagenta, kViolet, kBlue, kAzure, kCyan };
    if( layer == 0 ) return color_base[tile]+region;
    else return color_base[tile]-1-region;
    
  }
  
}

//_________________________________________________
void DrawRawDataClusterCharge_all()
{

  set_style(false);
  
  gStyle->SetOptStat(0);
  
  using run_pair_t = std::pair<int,int>;
  using run_map_t = std::vector<run_pair_t>;
//   run_map_t run_map =
//   {
//     {430, 21079},
//     {420, 21080},
//     {410, 21081},
//     {400, 21082},
//     {390, 21083},
//     {380, 21084},
//     {370, 21153},
//     {360, 21154},
//     {350, 21155},
//     {340, 21156},
//     {330, 21157},
//     {320, 21158}
//   };   
//   const TString pdffilename = "Figures/RawDataClusterCharge_phi-all.pdf";
  
//   run_map_t run_map =
//   {
//     {430, 21159},
//     {420, 21160},
//     {410, 21161},
//     {400, 21162},
//     {390, 21163},
//     {380, 21164},
//     {370, 21166},
//     {360, 21168},
//     {350, 21169},
//     {340, 21170},
//     {330, 21171},
//     {320, 21172}
//   };            
//   const TString pdffilename = "Figures/RawDataClusterCharge_z-all.pdf";

    
  // phi scan
  run_map_t run_map = 
  {
    { 430, 24341 },
    { 420, 24342 },
    { 410, 24343 },
    { 400, 24347 },
    { 390, 24348 },
    { 380, 24349 },
    { 370, 24350 },
    { 360, 24351 },
    { 350, 24352 },
    { 340, 24353 },
    { 330, 24354 },
    { 320, 24355 } 
  };
  const TString pdffilename = "Figures/RawDataClusterCharge_phi-all-new.pdf";

  PdfDocument pdfDocument( pdffilename );
  save_detector_names();

  // create TGraphErrors
  std::array<TGraphErrors*,64> tge_array = {};
  for( int ilayer = 0; ilayer <2; ++ilayer )
    for( int itile = 0; itile <8; ++itile )
  {
    const int detid = itile + 8*ilayer;
    const auto& detname = detector_names[detid];
    for( int iregion = 0; iregion<4; ++iregion )
    {
      int region_id = iregion +4*detid;
      const auto tgname = Form( "tge_%s_R%i", detname.c_str(), 4-iregion);
      const auto tgtitle =  Form( "cluster charge %s_R%i", detname.c_str(), 4-iregion);
      auto tge = new TGraphErrors();
      tge->SetName( tgname );
      tge->SetTitle( tgtitle );
      tge->GetXaxis()->SetTitle("Resist HV (V)");
      tge->GetYaxis()->SetTitle("#LTcluster charge#GT");
      tge->SetMarkerStyle( 20 );
      
      tge->SetMarkerColor( get_color( ilayer, itile, iregion ) );
      tge->SetLineColor( get_color( ilayer, itile, iregion ) );

      tge_array[region_id]=tge;
    }
  }
  
  // create tile TGraphErrors
  // create TGraphErrors
  std::array<TGraphErrors*,16> tge_tile_array = {};
  for( int ilayer = 0; ilayer <2; ++ilayer )
    for( int itile = 0; itile <8; ++itile )
  {
    const int detid = itile + 8*ilayer;
    const auto& detname = detector_names[detid];
    
    const auto tgname = Form( "tge_%s", detname.c_str());
    const auto tgtitle =  Form( "cluster charge %s", detname.c_str());
    auto tge = new TGraphErrors();
    tge->SetName( tgname );
    tge->SetTitle( tgtitle );
    tge->GetXaxis()->SetTitle("Resist HV (V)");
    tge->GetYaxis()->SetTitle("#LTcluster charge#GT");
    tge->SetMarkerStyle( 20 );
    
    tge->SetMarkerColor( get_color( ilayer, itile, 0 ) );
    tge->SetLineColor( get_color( ilayer, itile, 0 ) );
    
    tge_tile_array[detid]=tge;
  }
  
  // loop over runs
  int ipoint = 0;
  for( const auto& [hv,runnumber]:run_map )
  {
  
    const TString inputfilename = Form( "Rootfiles/RawDataClusters-%08i-0000.root", runnumber );
    
    FileManager fileManager( inputfilename );
    auto h_cluster_charge = static_cast<TH3*>( fileManager.GetHistogram( "h_cluster_charge" ) );
    
    for( int ilayer = 0; ilayer <2; ++ilayer )
      for( int itile = 0; itile <8; ++itile )
    {
      const int detid = itile + 8*ilayer;
      for( int iregion = 0; iregion<4; ++iregion )
      {
        int region_id = iregion +4*detid;
        h_cluster_charge->GetXaxis()->SetRange(region_id+1, region_id+1);
        
        // only select clusters of size 2 or more
        // h_cluster_charge->GetYaxis()->SetRange( 3, h_cluster_charge->GetNbinsY()+1 );
        
        // project
        auto h = std::unique_ptr<TH1>( static_cast<TH1*>( h_cluster_charge->Project3D( "z" ) ) );
        const auto hname = Form( "h_%s_R%i", detector_names[detid].c_str(), 4-iregion );
        const auto htitle = Form( "cluster charge %s_R%i", detector_names[detid].c_str(), 4-iregion );
        h->SetName( hname );
        h->SetTitle( htitle );
      
        const auto mean = h->GetMean();
        const auto error = h->GetMeanError();
        
        tge_array[region_id]->SetPoint(ipoint, hv, mean );
        tge_array[region_id]->SetPointError( ipoint, 0, error );
        
      }
      
      {
        // project over entire detector
        h_cluster_charge->GetXaxis()->SetRange(4*detid+1, 4*detid+4);
        auto h = std::unique_ptr<TH1>( static_cast<TH1*>( h_cluster_charge->Project3D( "z" ) ) );
        const auto hname = Form( "h_%s", detector_names[detid].c_str() );
        const auto htitle = Form( "cluster charge %s", detector_names[detid].c_str() );
        h->SetName( hname );
        h->SetTitle( htitle );
        
        const auto mean = h->GetMean();
        const auto error = h->GetMeanError();
        
        tge_tile_array[detid]->SetPoint(ipoint, hv, mean );
        tge_tile_array[detid]->SetPointError( ipoint, 0, error );
      }
      
    }
    
    ++ipoint;
    
  }
  
  if( true )
  {
    auto cv = new TCanvas( "cv_all", "cv_all", 1200, 1200 );
    cv->Divide( 4, 4, 0.002, 0.002 );
    
    // one canvas per region
    for( int ilayer = 0; ilayer <2; ++ilayer )
      for( int itile = 0; itile <8; ++itile )
    {
      const int detid = itile + 8*ilayer;
      const auto& detname = detector_names[detid];

      cv->cd( detid+1 );
      
      const auto hname = Form( "h_%i", detid );
      auto h = new TH1I( hname, "", 100, 310, 440 );
      h->SetMinimum( 90 );
      h->SetMaximum( 1200 );
      h->GetXaxis()->SetTitle( "Resist HV (V)" );
      h->GetYaxis()->SetTitle( "cluster charge" );
      h->Draw();
      gPad->SetLogy( true );
      
      // legend
      auto legend = new TLegend( 0.15, 0.6, 0.4, 0.85, "", "NDC" );
      legend->SetFillColor(0);
      legend->SetFillStyle(0);
      legend->SetBorderSize(0);
      legend->Draw();
      
      for( int iregion = 0; iregion<4; ++iregion )
      {
        int region_id = iregion +4*detid;
        tge_array[region_id]->Draw("P");
        legend->AddEntry( tge_array[region_id], Form( "%s_R%i", detname.c_str(), 4-iregion ), "PL" );
      }
    }
    
    pdfDocument.Add( cv );
    
  }
 
  if( true )
  {
    auto cv = new TCanvas( "cv_tile_all", "cv_all", 1200, 1200 );
    cv->Divide( 4, 4, 0.002, 0.002 );
    
    // one canvas per region
    for( int ilayer = 0; ilayer <2; ++ilayer )
      for( int itile = 0; itile <8; ++itile )
    {
      const int detid = itile + 8*ilayer;
      const auto& detname = detector_names[detid];

      cv->cd( detid+1 );
      
      const auto hname = Form( "h_tile_%i", detid );
      auto h = new TH1I( hname, "", 100, 310, 440 );
      h->SetMinimum( 90 );
      h->SetMaximum( 1200 );
      h->GetXaxis()->SetTitle( "Resist HV (V)" );
      h->GetYaxis()->SetTitle( "cluster charge" );
      h->Draw();
      gPad->SetLogy( true );
      
      tge_tile_array[detid]->Draw("P");
      Draw::PutText( 0.15, 0.8, detname.c_str(), 0.1);
    }
    
    pdfDocument.Add( cv );
    
  }
}
