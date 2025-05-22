#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

#include <TH1.h>
#include <TF1.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>

//_____________________________________________________________________________
void DrawRawDataTiming_driftscan()
{
  MicromegasMapping mapping;
  const TString pdfFile = "Figures/DrawRawDataTiming_driftscan-all-v3.pdf";
  PdfDocument pdfDocument( pdfFile );

  
  // map HV to run number
  using runmap_t = std::map<int,int>;
//   runmap_t run_map = 
//   {
//     { 50, 7409 },
//     { 100, 7410 },
//     { 150, 7411 },
//     { 200, 7412 },
//     { 250, 7413 },
//     { 300, 7414 }
//   };
    
//   runmap_t run_map = 
//   {
//     { 50, 7450 },
//     { 100, 7451 },
//     { 150, 7449 },
//     { 200, 7452 },
//     { 250, 7453 },
//     { 300, 7454 },
//     { 350, 7455 },
//     { 400, 7456 }
//   };

  runmap_t run_map = 
  {
    { 50,  9429 },
    { 100, 9430 },
    { 150, 9431 },
    { 200, 9432 },
    { 250, 9433 },
    { 300, 9434 },
    { 350, 9435 },
    { 400, 9436 }
  };

  // open all TFiles
  using tfilemap_t = std::map<int,TFile*>;
  tfilemap_t tfile_map;
  
  for( const auto& [hv, runnumber]:run_map )
  {
    const auto filename = Form( "Rootfiles/RawDataTiming-%08i-0000.root", runnumber );
    auto tfile = TFile::Open( filename );
    tfile_map.insert( std::make_pair( hv, tfile ) );
  }
  
//   int ilayer = 0;
//   int itile = 1;
//   int iregion = 0;
//   int itile = 1;
  for( int ilayer = 0; ilayer < 2; ++ilayer )
    for( int itile = 0; itile < 8; ++itile )
    for( int iregion = 0; iregion < 4; ++iregion )
  {
    const auto cvname = Form( "cv_%i_%i_%i", ilayer, itile, iregion );
    TCanvas* cv = new TCanvas( cvname, cvname, 1200, 900 );
    cv->Divide( 4, 3 );
    
    // get detector name
    const int layer = 55+ilayer;
    const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
    const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
    const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
    
    int ipad = 0;
    for( const auto& [hv, tfile]: tfile_map )
    {
      
      const auto runnumber = run_map[hv];
      
      
      // create histogram name
      const auto hname = Form( "h_%i_%i_%i", ilayer, itile, iregion );
      auto h = static_cast<TH1*>( tfile->Get(hname) );
      if( !h ) 
      {
        std::cout << "hv: " << hv << " histogram " << hname << " not found" << std::endl;
        continue;
      }
      
      std::cout << "HV: " << hv << " Detector name: " << name << " histogram: " << h->GetTitle() << std::endl;
      
      // change histogram name and title
      const auto hname_new = Form( "%s_%iV", h->GetName(), hv );
      const auto htitle_new = Form( "%s (%iV) (%i)", h->GetTitle(), hv, runnumber );
      h->SetName( hname_new );
      h->SetTitle( htitle_new );
      cv->cd( ++ipad );
      h->Draw("colz" );
      
    }
    pdfDocument.Add( cv );
    
  }

  
  
}
    
