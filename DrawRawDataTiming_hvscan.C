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
void DrawRawDataTiming_hvscan()
{
  MicromegasMapping mapping;
  // const TString pdfFile = "Figures/DrawRawDataTiming.pdf";
  const TString pdfFile = "./Figures/DrawRawDataTiming_hvscan.pdf";
  PdfDocument pdfDocument( pdfFile );

  
  // map HV to run number
// /bbox/commissioning/TPOT/physics/TPOT_ebdc39_physics-00007362-0000.prdf All nominal HV 1 minute
// /bbox/commissioning/TPOT/physics/TPOT_ebdc39_physics-00007363-0000.prdf All HV at 370 1 minute
// /bbox/commissioning/TPOT/physics/TPOT_ebdc39_physics-00007364-0000.prdf All HV at 350
// /bbox/commissioning/TPOT/physics/TPOT_ebdc39_physics-00007365-0000.prdf All HV at 320 (= Safe)
// /bbox/commissioning/TPOT/physics/TPOT_ebdc39_physics-00007366-0000.prdf All HV at 330
// /bbox/commissioning/TPOT/physics/TPOT_ebdc39_physics-00007367-0000.prdf All HV at 340
// /bbox/commissioning/TPOT/physics/TPOT_ebdc39_physics-00007368-0000.prdf All HV at 360
// /bbox/commissioning/TPOT/physics/TPOT_ebdc39_physics-00007369-0000.prdf Nominal
// /bbox/commissioning/TPOT/physics/TPOT_ebdc39_physics-00007370-0000.prdf 2 minutes 420 V max
// /bbox/commissioning/TPOT/physics/TPOT_ebdc39_physics-00007371-0000.prdf 2 minutes 410 V max
// /bbox/commissioning/TPOT/physics/TPOT_ebdc39_physics-00007372-0000.prdf 2 minutes 400 V max
// /bbox/commissioning/TPOT/physics/TPOT_ebdc39_physics-00007375-0000.prdf 2 minutes 390 V max
// /bbox/commissioning/TPOT/physics/TPOT_ebdc39_physics-00007376-0000.prdf 2 minutes 380 V max
  using runmap_t = std::map<int,int>;
  runmap_t run_map = {
    { 430, 7362 },
//     { 430, 7369 },
    { 370, 7363 },
    { 350, 7364 },
    { 320, 7365 },
    { 330, 7366 },
    { 340, 7367 },
    { 360, 7368 },
    { 420, 7370 },
    { 410, 7371 },
    { 400, 7372 },
    { 390, 7375 },
    { 380, 7376 } 
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
    
