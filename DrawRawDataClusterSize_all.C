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

}

//_________________________________________________
void DrawRawDataClusterSize_all()
{

  using run_pair_t = std::pair<int,int>;
  using run_map_t = std::vector<run_pair_t>;
  run_map_t run_map =
  {
    {430, 9439},
    {420, 9440},
    {410, 9441},
//     {410, 9442},
    {400, 9443},
    {390, 9444},
    {380, 9446},
    {370, 9447},
    {360, 9448},
    {350, 9449},
    {340, 9450},
    {330, 9451},
    {320, 9452}
  };

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
      const auto tgtitle =  Form( "cluster size %s_R%i", detname.c_str(), 4-iregion);
      auto tge = new TGraphErrors();
      tge->SetName( tgname );
      tge->SetTitle( tgtitle );
      tge->GetXaxis()->SetTitle("HV (V)");
      tge->GetYaxis()->SetTitle("#LTcluster size#GT");
      tge->SetMarkerStyle( 20 );
      tge->SetMarkerColor(1);
      
      tge_array[region_id]=tge;
    }
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
        auto h = std::unique_ptr<TH1>( static_cast<TH1*>( h_cluster_charge->Project3D( "y" ) ) );
        const auto hname = Form( "h_%s_R%i", detector_names[detid].c_str(), 4-iregion );
        const auto htitle = Form( "cluster size %s_R%i", detector_names[detid].c_str(), 4-iregion );
        h->SetName( hname );
        h->SetTitle( htitle );
      
        const auto mean = h->GetMean();
        const auto error = h->GetMeanError();
        
        tge_array[region_id]->SetPoint(ipoint, hv, mean );
        tge_array[region_id]->SetPointError( ipoint, 0, error );
        
      }
    }
    
    ++ipoint;
    
  }
  
  // make plot
  const TString pdffilename = "Figures/RawDataClusterSize-all.pdf";
  PdfDocument pdfDocument( pdffilename );

  for( int ilayer = 0; ilayer <2; ++ilayer )
    for( int itile = 0; itile <8; ++itile )
  {
    const int detid = itile + 8*ilayer;
    auto cvname = Form( "cv_%i", detid);
    auto cv = std::make_unique<TCanvas>( cvname, cvname, 900,500 );
    cv->Divide( 2, 2 );
    for( int iregion = 0; iregion<4; ++iregion )
    {
      cv->cd( iregion+1);
      int region_id = iregion +4*detid;
      tge_array[region_id]->Draw("AP");
    }
    pdfDocument.Add( cv.get() );
  }
  
}
