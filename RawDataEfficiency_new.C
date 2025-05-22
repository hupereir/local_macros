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
eff_vector_t get_detector_efficiencies( const TString& path, int runnumber )
{
  eff_vector_t out(16);

  const TString inputFile = Form( "%s/dst_eval-%08i-0000-full.root", path.Data(), runnumber );
  // const TString inputFile = Form( "%s/dst_eval-%08i-0000.root", path.Data(), runnumber );
  std::cout << "RawDataEfficiency_single - inputFile: " << inputFile << std::endl;

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  for( int ilayer = 0; ilayer < 2; ++ilayer )
    for( int itile = 0; itile < 8; ++itile )
  {

    if( itile != 3 ) continue;


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
      << " ref: " << h->GetEntries()
      << " found: " << h->Integral( 2, h->GetNbinsX()+1)
      << " efficiency: " << eff << "+/-" << err
      << std::endl;

    out[detid] = {eff,err};
    delete h;
  }
  return out;
}

//__________________________________________________________________________________
void RawDataEfficiency_new()
{
  save_detector_names();

  // phi scan (D400)
  using run_map_t = std::vector<std::pair<int,int>>;
//   const TString path = "DST/CONDOR_CombinedDataEvaluation";
//   const TString suffix = "";
//   const TString tag( "_phi_2024_06_27-D400" );
//   run_map_t run_map =
//   {
//     { 410, 46810 },
//     { 320, 46811 },
//     { 330, 46812 },
//     { 340, 46813 },
//     { 350, 46814 },
//     { 360, 46815 },
//     { 370, 46818 },
//     { 380, 46819 },
//     { 380, 46821 },
//     { 390, 46822 },
//     { 400, 46823 },
//     { 410, 46824 },
//     { 420, 46825 },
//     { 430, 46826 },
//     { 440, 46827 },
//     { 450, 46828 },
//     {0,0}
//   };

//   const TString path = "DST/CONDOR_CombinedDataEvaluation";
//   const TString suffix = "";
//   const TString tag( "_z_2024_06_27-D400" );
//   run_map_t run_map =
//   {
//     { 410, 46829 },
//     { 320, 46831 },
//     { 330, 46835 },
//     { 340, 46836 },
//     {0,0}
//   };

//   const TString path = "DST/CONDOR_CombinedDataEvaluation_new_calib_zs";
//   const TString suffix = "_new_calib_zs";
//
//   const TString path = "DST/CONDOR_CombinedDataEvaluation";
//   const TString suffix = "";
//
//   const TString tag( "_phi_2024_07_08-D400" );
//   run_map_t run_map =
//   {
//     { 410, 47583 },
//     { 320, 47584 },
//     { 330, 47585 },
//     { 340, 47586 },
//     { 350, 47587 },
//     { 360, 47589 },
//     { 370, 47590 },
//     { 380, 47591 },
//     { 390, 47594 },
//     { 400, 47595 },
//     { 420, 47602 },
//     { 430, 47603 },
//     { 440, 47604 },
//     { 450, 47605 },
//     { 410, 47607 },
//     {0,0}
//   };
//
//   const TString path = "DST/CONDOR_CombinedDataEvaluation_new_calib_zs";
//   const TString suffix = "_new_calib_zs";
// //   const TString path = "DST/CONDOR_CombinedDataEvaluation";
// //   const TString suffix = "";
//   const TString tag( "_z_2024_07_08-D400" );
//   run_map_t run_map =
//   {
//     { 410, 47607 },
//     { 430, 47608 },
//     { 420, 47620 },
//     { 410, 47621 },
//     { 400, 47622 },
//     { 390, 47624 },
//     { 380, 47625 },
//     { 370, 47626 },
//     { 360, 47629 },
//     { 350, 47630 },
//     { 340, 47631 },
//     { 330, 47632 },
//     { 320, 47633 },
//     {0,0}
//   };


// //   const TString path = "DST/CONDOR_CombinedDataEvaluation";
// //   const TString suffix = "";
//   const TString suffix = "_new_calib_zs";
//   const TString tag( "_phi_2024_08_09-D400" );
//   run_map_t run_map =
//   {
//     { 410, 50963 },
//     { 450, 50966 },
//     { 440, 50968 },
//     { 440, 50970 },
//     { 430, 50971 },
//     { 420, 50972 },
//     { 410, 50974 },
//     { 400, 50975 },
//     { 390, 50976 },
//     { 380, 50977 },
//     { 370, 50978 },
//     { 360, 50980 },
//     { 360, 50982 },
//     { 350, 50983 },
//     { 340, 50984 },
//     { 330, 50985 },
//     { 320, 50986 },
//     {0,0}
//   };


  const TString path = "DST/CONDOR_CombinedDataEvaluation";
  // const TString path = "DST/";
  const TString suffix = "_new_calib_zs";
  const TString tag( "_2024_08_28-D400" );
  run_map_t run_map =
  {
    { 410, 50963 },
    {0,0}
  };
  const TString rootfilename = Form("Rootfiles/RawDataEfficiency_new%s%s.root",tag.Data(), suffix.Data());
  const TString pdffilename = Form("Figures/RawDataEfficiency_new%s%s.pdf",tag.Data(), suffix.Data());

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

    // store efficiency and error in array
    std::vector<double>eff_array;
    std::vector<double>err_array;

    if( !hv ) continue;
    const auto efficiencies = get_detector_efficiencies( path, runnumber );
    for( int detid = 0; detid<16; ++detid )
    {
      const auto& [eff,err] = efficiencies[detid];
      eff_array.push_back(eff);
      err_array.push_back(err);

      tge_tile_array[detid]->SetPoint( ipoint, hv, eff );
      tge_tile_array[detid]->SetPointError( ipoint, 0, err );
    }
    ++ipoint;

    // print
    Stream::PrintVector( "eff", eff_array );
    Stream::PrintVector( "err", err_array );

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
      h->SetStats(0);
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
