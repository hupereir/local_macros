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

//__________________________________________________________________________________
using eff_pair_t = std::pair<double,double>;
using eff_vector_t = std::vector<eff_pair_t>;
eff_vector_t get_detector_efficiencies( int runnumber )
{
  eff_vector_t out(16);

  const TString inputFile = Form( "DST/CONDOR_CombinedDataEvaluation/dst_eval-%08i-0000-full.root", runnumber );
  std::cout << "RawDataEfficiency_single - inputFile: " << inputFile << std::endl;

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  for( int ilayer = 0; ilayer < 2; ++ilayer )
    for( int itile = 0; itile < 8; ++itile )
  {
    const int detid = itile + 8*ilayer;
    const int detid_ref = (ilayer==0) ? detid+8:detid-8;

    const auto var = Form( "n_detector_clusters[%i]", detid );
    const TCut cut = Form( "n_detector_clusters[%i]==1", detid_ref );
    const TCut min_cluster_charge_cut = Form( "min_cluster_charge[%i]>200", detid_ref );
    // const TCut min_cluster_charge_cut = Form( "min_cluster_charge[%i]>800", detid_ref );

    const auto hname = Form( "h_%i", detid );
    auto h = new TH1I( hname, "", 50, 0, 50 );
    Utils::TreeToHisto( tree, hname, var, cut && min_cluster_charge_cut, false );

    const double eff = h->Integral( 2, h->GetNbinsX()+1)/h->GetEntries();
    const double err = std::sqrt( eff*(1.0-eff)/h->GetEntries() );
    std::cout
      << "get_detector_efficiencies -"
      << " runnumber: " << runnumber
      << " detid: " << detid
      << " name: " << detector_names[detid]
      << " efficiency: " << eff << "+/-" << err
      << std::endl;

    out[detid] = {eff,err};
    delete h;
  }
  return out;
}

//__________________________________________________________________________________
void RawDataEfficiency_for_bade()
{
  save_detector_names();

  // phi scan (D400)
  using run_map_t = std::vector<std::pair<int,int>>;
  const TString tag( "_phi_2023-D400" );
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

  const TString rootfilename = Form("Rootfiles/RawDataEfficiency_new%s.root",tag.Data());
  const TString pdffilename = Form("Figures/RawDataEfficiency_new%s.pdf",tag.Data());

  std::cout << "RawDataEfficiency - rootfilename: " << rootfilename << std::endl;
  std::cout << "RawDataEfficiency - pdffilename: " << pdffilename << std::endl;

  PdfDocument pdfdocument( pdffilename );
  RootFile rootfile( rootfilename );

  // create TGraphErrors
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
  }

  int ipoint = 0;
  for( const auto& [hv,runnumber]:run_map )
  {
    if( !hv ) continue;
    const auto efficiencies = get_detector_efficiencies( runnumber );
    for( int detid = 0; detid<16; ++detid )
    {
      const auto& [eff,err] = efficiencies[detid];
      tge_tile_array[detid]->SetPoint( ipoint, hv, eff );
      tge_tile_array[detid]->SetPointError( ipoint, 0, err );
    }
    ++ipoint;
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
      auto h = new TH1I( hname, "", 100, 310, 460 );
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
