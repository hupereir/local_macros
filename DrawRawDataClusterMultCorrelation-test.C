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
    for( int itile = 0; itile < 8; ++itile )
    {
      const int layer = 55;
      const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
      detector_names.push_back( std::move( name.substr(0,3) ) );
    }
  }

}

//_________________________________________________
void DrawRawDataClusterMultCorrelation_test()
{
  int runNumber = 9443;  

  const TString inputfilename = Form( "Rootfiles/RawDataClusters-%08i-0000.root", runNumber );
  const TString pdffilename = Form( "Figures/RawDataClusterMultCorrelation-%08i-0000_test.pdf", runNumber );

  // pdffile
  PdfDocument pdfDocument( pdffilename );

  FileManager fileManager( inputfilename );
  auto h3d = static_cast<TH3*>( fileManager.GetHistogram( "h_cluster_mult_correlation" ) );
  
  save_detector_names();
  
  if( true )
  {
    for( int itile = 0; itile <8; ++itile )
    {
      const int detid = itile;
      auto cvname = Form( "cv_%i", detid);
      auto cv = std::make_unique<TCanvas>( cvname, cvname, 900,500 );
      
      h3d->GetXaxis()->SetRange(detid+1, detid+1);
      int min_cluster = 1;
      int max_cluster = 10;
//       h3d->GetYaxis()->SetRangeUser(1, 50 );
//       h3d->GetZaxis()->SetRangeUser(1, 50 );

      auto h2d = static_cast<TH2*>( h3d->Project3D( "zy" ) );
      {
        const auto hname = Form( "h2d_%s", detector_names[detid].c_str() );
        const auto htitle = Form( "cluster multiplicity %s", detector_names[detid].c_str() );
        h2d->SetName( hname );
        h2d->SetTitle( htitle );
      }
      
      {
        const auto hname = Form( "h1d_%s", detector_names[detid].c_str() );
        const auto htitle = Form( "cluster mult (#phi-z) %s", detector_names[detid].c_str() );
        auto h1d = new TH1I( hname, htitle, 2*max_cluster, -max_cluster, max_cluster );
        h1d->GetXaxis()->SetTitle( "cluster mult (#phi-z)" );
        
        for( int ix = 1; ix < h2d->GetNbinsX(); ++ix )
          for( int iy = 1; iy < h2d->GetNbinsY(); ++iy )
        {
          double mx = h2d->GetXaxis()->GetBinCenter( ix+1); 
          double my = h2d->GetYaxis()->GetBinCenter( iy+1); 
          double content = h2d->GetBinContent( ix+1, iy+1 );
          h1d->Fill( mx-my, content );
        }
        
        h1d->Draw( "hist" );
        pdfDocument.Add( cv.get() );
        
      }
    }
  }

}
