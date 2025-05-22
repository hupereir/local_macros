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

//___________________________________________________________________-
void RawDataEfficiency_single( int runnumber = 47604)
{
  save_detector_names();

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
      << "RawDataEfficiency -"
      << " runnumber: " << runnumber
      << " detid: " << detid
      << " name: " << detector_names[detid]
      << " efficiency: " << eff << "+/-" << err
      << std::endl;
  }
}
