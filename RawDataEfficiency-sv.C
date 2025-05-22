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

enum cut_id_t
{
  no_cut,
  cluster_charge_cut,
  cluster_size_cut,
  cluster_charge_size_cut,
  sigma4_no_cut,
  sigma4_cluster_charge_cut,
  sigma4_cluster_size_cut,
  sigma4_cluster_charge_size_cut
};

void RawDataEfficiency( cut_id_t cut_id = sigma4_cluster_charge_size_cut )
{
  
  std::cout << "RawDataEfficiency - cut_id: " << cut_id << std::endl;

  TString suffix = "";
  switch( cut_id )
  {
    case no_cut: suffix = "_nocut"; break;
    case cluster_charge_cut: suffix = "_charge_cut"; break;
    case cluster_size_cut: suffix = "_size_cut"; break;
    case cluster_charge_size_cut: suffix = "_charge_size_cut"; break;
    case sigma4_no_cut: suffix = "_4sigma"; break;
    case sigma4_cluster_charge_cut: suffix = "_4sigma_charge_cut"; break;
    case sigma4_cluster_size_cut: suffix = "_4sigma_size_cut"; break;
    case sigma4_cluster_charge_size_cut: suffix = "_4sigma_charge_size_cut"; break;
  }
    
    
  using run_map_t = std::vector<std::pair<int,int>>;
  save_detector_names();
  
//   // phi scan (D300)
//   run_map_t run_map = 
//   {
//     { 430, 21079 },
//     { 420, 21080 },
//     { 410, 21081 },
//     { 400, 21082 },
//     { 390, 21083 },
//     { 380, 21084 },
//     { 370, 21153 },
//     { 360, 21154 },
//     { 350, 21155 },
//     { 340, 21156 },
//     { 330, 21157 },
//     { 320, 21158 } 
//   };
//     
//   const TString rootfilename = "Rootfiles/RawDataEfficiency_phi-D300.root";
//   const TString pdffilename = "Figures/RawDataEfficiency_phi-D300.pdf";

//   // z scan (D300)
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
//   
//   const TString rootfilename = "Rootfiles/RawDataEfficiency_z-D300.root";
//   const TString pdffilename = "Figures/RawDataEfficiency_z-D300.pdf";
// 
//   // phi scan (D400)
//   run_map_t run_map = 
//   {
//     { 430, 24341 },
//     { 420, 24342 },
//     { 410, 24343 },
//     { 400, 24347 },
//     { 390, 24348 },
//     { 380, 24349 },
//     { 370, 24350 },
//     { 360, 24351 },
//     { 350, 24352 },
//     { 340, 24353 },
//     { 330, 24354 },
//     { 320, 24355 } 
//   };
//     
//   const TString rootfilename = Form( "Rootfiles/RawDataEfficiency_phi-D400%s.root", suffix.Data() );
//   const TString pdffilename = Form( "Figures/RawDataEfficiency_phi-D400%s.pdf", suffix.Data() );
  
  // z scan (D400)
  run_map_t run_map =
  {
    {430, 24358},
    {420, 24359},
    {410, 24360},
    {400, 24365},
    {390, 24366},
    {380, 24367},
    {370, 24368},
    {360, 24369},
    {350, 24370},
    {340, 24371},
    {330, 24372},
    {320, 24373}
  };            
  
  const TString rootfilename = Form( "Rootfiles/RawDataEfficiency_z-D400%s.root", suffix.Data() );
  const TString pdffilename = Form( "Figures/RawDataEfficiency_z-D400%s.pdf", suffix.Data() );
  
    
  std::cout << "RawDataEfficiency - rootfilename: " << rootfilename << std::endl;
  std::cout << "RawDataEfficiency - pdffilename: " << pdffilename << std::endl;

  PdfDocument pdfdocument( pdffilename );
  RootFile rootfile( rootfilename ); 

  // create TGraphErrors
  std::array<TGraphErrors*,64> tge_region_array = {};
  std::array<TGraphErrors*,16> tge_tile_array = {};
  for( int ilayer = 0; ilayer <2; ++ilayer )
    for( int itile = 0; itile <8; ++itile )
  {
    const int detid = itile + 8*ilayer;
    const auto& detname = detector_names[detid];
    
    const auto tgname = Form( "tge_%s", detname.c_str());
    const auto tgtitle =  Form( "Efficiency %s", detname.c_str());
    auto tge = new TGraphErrors();
    tge->SetName( tgname );
    tge->SetTitle( tgtitle );
    tge->GetXaxis()->SetTitle("Resist HV (V)");
    tge->GetYaxis()->SetTitle("efficiency");
    tge->SetMarkerStyle( 20 );
    
    tge->SetMarkerColor( get_color( ilayer, itile, 0 ) );
    tge->SetLineColor( get_color( ilayer, itile, 0 ) );
    
    tge_tile_array[detid]=tge;
    rootfile.Add(tge);
  
    for( int iregion = 0; iregion < 4; ++iregion )
    {
      const int regionid = iregion +4*detid;
      const auto regionname  = Form( "%s_R%i", detname.c_str(), 4-iregion );
      const auto tgname = Form( "tge_%s", regionname);
      const auto tgtitle =  Form( "Efficiency %s", regionname);
      auto tge = new TGraphErrors();
      tge->SetName( tgname );
      tge->SetTitle( tgtitle );
      tge->GetXaxis()->SetTitle("Resist HV (V)");
      tge->GetYaxis()->SetTitle("efficiency");
      tge->SetMarkerStyle( 20 );
      
      tge->SetMarkerColor( get_color( ilayer, itile, iregion ) );
      tge->SetLineColor( get_color( ilayer, itile, iregion ) );
      
      tge_region_array[regionid]=tge;
      rootfile.Add(tge);
    }
    
  }
  
  int ipoint = 0;
  for( const auto& [hv,runnumber]:run_map )
  {
    const bool is_4sigma = 
      cut_id == sigma4_no_cut || cut_id == sigma4_cluster_charge_cut || 
      cut_id == sigma4_cluster_size_cut || cut_id == sigma4_cluster_charge_size_cut;
    const TString rootfilename = is_4sigma ?  Form( "Rootfiles/RawDataClusterTree-%08i-0000-full_4sigma.root", runnumber ) : Form( "Rootfiles/RawDataClusterTree-%08i-0000-full.root", runnumber );
    std::cout << "RawDataEfficiency - rootfilename: " << rootfilename << std::endl;
    auto tfile = std::unique_ptr<TFile>( TFile::Open( rootfilename, "READ" ) );
    auto tree = static_cast<TTree*>(tfile->Get("T"));
        
    for( int ilayer = 0; ilayer < 2; ++ilayer )
      for( int itile = 0; itile < 8; ++itile )
    {
      const int detid = itile + 8*ilayer;
      const int detid_ref = (ilayer==0) ? detid+8:detid-8;
      
      const auto var = Form( "n_detector_clusters[%i]", detid );
      TCut cut = Form( "n_detector_clusters[%i]==1", detid_ref );
      /*
       * for detid 5 (NEIP), one hast to also require that NEIZ_R3 does not have any clusters:
       * this region is indeed noisy, which systematically reduce the efficiency measured on NEIP
       */
      if( detid == 5 ) cut = cut && Form( "n_region_clusters[%i]==0", 4*detid_ref+1 );
      
      TCut min_cluster_charge_cut = Form( "min_cluster_charge[%i]>200", detid_ref ); // that looks pretty ok
      TCut min_cluster_size_cut = Form( "min_cluster_size[%i]>=2", detid_ref );
      switch( cut_id )
      {
        default: break;
        case cluster_charge_cut:
        case sigma4_cluster_charge_cut:
        cut = cut&&min_cluster_charge_cut;
        break;
        
        case cluster_size_cut:
        case sigma4_cluster_size_cut:
        cut = cut&&min_cluster_size_cut;
        break;
        
        case cluster_charge_size_cut:
        case sigma4_cluster_charge_size_cut:
        cut = cut&&min_cluster_size_cut&&min_cluster_charge_cut;
        break;
      }
      
      const auto hname = Form( "h_%i_%i", detid, hv );
      auto h = new TH1I( hname, "", 50, 0, 50 );
      
      Utils::TreeToHisto( tree, hname, var, cut, false );
    
      const double eff = 1. - double( h->GetBinContent(1) )/h->GetEntries();
      const double err = std::sqrt( eff*(1.0-eff)/h->GetEntries() );
      
      const double eff2 = h->Integral( 2, h->GetNbinsX()+1)/h->GetEntries();
      std::cout << "detid -eff: " << eff << " eff2: " << eff2 << std::endl;
      
      tge_tile_array[detid]->SetPoint( ipoint, hv, eff );
      tge_tile_array[detid]->SetPointError( ipoint, 0, err );
    
      delete h;
      
      for( int iregion =0; iregion < 4; ++iregion )
      {
        const int regionid = iregion+4*detid;
        const auto region_var = Form( "n_region_clusters[%i]", regionid );
        
        const auto hname = Form( "h_r%i_%i", regionid, hv );
        auto h = new TH1I( hname, "", 50, 0, 50 );
        
        Utils::TreeToHisto( tree, hname, region_var, cut, false );
    
        const double eff = 1. - double( h->GetBinContent(1) )/h->GetEntries();
        const double err = std::sqrt( eff*(1.0-eff)/h->GetEntries() );
      
        tge_region_array[regionid]->SetPoint( ipoint, hv, eff );
        tge_region_array[regionid]->SetPointError( ipoint, 0, err );
        
        delete h;
        
      }
    
    }
    ++ipoint;
  }
  
  // make plot
  if( true )
  {
    auto cv = new TCanvas( "cv_region_all", "cv_all", 1200, 1200 );
    cv->Divide( 4, 4, 0.002, 0.002 );
    
    // one canvas per region
    for( int ilayer = 0; ilayer <2; ++ilayer )
      for( int itile = 0; itile <8; ++itile )
    {
      const int detid = itile + 8*ilayer;
      const auto& detname = detector_names[detid];

      cv->cd( detid+1 );
      
      const auto hname = Form( "h_region_%i", detid );
      auto h = new TH1I( hname, "", 100, 310, 440 );
      h->SetMinimum( 0 );
      h->SetMaximum( 0.5 );
      h->GetXaxis()->SetTitle( "Resist HV (V)" );
      h->GetYaxis()->SetTitle( "efficiency" );
      h->Draw();
    
      // legend
      auto legend = new TLegend( 0.15, 0.6, 0.4, 0.85, "", "NDC" );
      legend->SetFillColor(0);
      legend->SetFillStyle(0);
      legend->SetBorderSize(0);
      legend->Draw();
      
      for( int iregion = 0; iregion < 4; ++iregion )
      {
        int regionid = iregion + 4*detid; 
        tge_region_array[regionid]->Draw("P");
        legend->AddEntry( tge_region_array[regionid], Form( "%s_R%i", detname.c_str(), 4-iregion ), "PL" );
      }
    
    
    }
    pdfdocument.Add( cv );
    
  }
  
  // make plot
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
      h->SetMinimum( 0 );
      h->SetMaximum( 1. );
      h->GetXaxis()->SetTitle( "Resist HV (V)" );
      h->GetYaxis()->SetTitle( "efficiency" );
      h->Draw();
      
      tge_tile_array[detid]->Draw("P");
      Draw::PutText( 0.15, 0.8, detname.c_str(), 0.1);
    }
    
    pdfdocument.Add( cv );
    
  }

}

void process_all()
{
  for( const cut_id_t& cut_id: {
    no_cut,
    cluster_charge_cut,
    cluster_size_cut,
    cluster_charge_size_cut,
    sigma4_no_cut,
    sigma4_cluster_charge_cut,
    sigma4_cluster_size_cut,
    sigma4_cluster_charge_size_cut
  } ) { 
    
    std::cout << std::endl << std::endl << std::endl;
    RawDataEfficiency( cut_id ); 

  }

}
