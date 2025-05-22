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
void DrawRawDataClusterSize()
{
  
  gStyle->SetOptStat(111111);
  
  // int runNumber = 9439;  
  int runNumber = 9440;  
  // int runNumber = 13010;  

  const TString inputfilename = Form( "Rootfiles/RawDataClusters-%08i-0000.root", runNumber );
  const TString pdffilename = Form( "Figures/RawDataClusterSize-%08i-0000.pdf", runNumber );

  // pdffile
  PdfDocument pdfDocument( pdffilename );

  FileManager fileManager( inputfilename );
  auto h_cluster_charge = static_cast<TH3*>( fileManager.GetHistogram( "h_cluster_charge" ) );

  save_detector_names();
  
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
      h_cluster_charge->GetXaxis()->SetRange(region_id+1, region_id+1);
      auto h = static_cast<TH1*>( h_cluster_charge->Project3D( "y" ) );
      
      const auto hname = Form( "h_%s_R%i", detector_names[detid].c_str(), 4-iregion );
      const auto htitle = Form( "cluster size %s_R%i", detector_names[detid].c_str(), 4-iregion );
      h->SetName( hname );
      h->SetTitle( htitle );
      
      h->Draw();
      gPad->SetLogz(true);
    }
    pdfDocument.Add( cv.get() );
  }
}
