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
  cluster_strip_cut,
  cluster_charge_cut,
  cluster_size_cut,
  cluster_charge_size_cut,
  sigma4_no_cut,
  sigma4_cluster_charge_cut,
  sigma4_cluster_size_cut,
  sigma4_cluster_charge_size_cut
};

void RawDataEfficiency_2d( int runnumber = 24341, cut_id_t cut_id = no_cut )
{

  save_detector_names();

  TString suffix = "";
  switch( cut_id )
  {
    case no_cut: suffix = "_nocut"; break;
    case cluster_strip_cut: suffix = "_strip_cut"; break;
    case cluster_charge_cut: suffix = "_charge_cut"; break;
    case cluster_size_cut: suffix = "_size_cut"; break;
    case cluster_charge_size_cut: suffix = "_charge_size_cut"; break;

    case sigma4_no_cut: suffix = "_4sigma"; break;
    case sigma4_cluster_charge_cut: suffix = "_4sigma_charge_cut"; break;
    case sigma4_cluster_size_cut: suffix = "_4sigma_size_cut"; break;
    case sigma4_cluster_charge_size_cut: suffix = "_4sigma_charge_size_cut"; break;

//     case sigma4_no_cut: suffix = "_4.5sigma"; break;
//     case sigma4_cluster_charge_cut: suffix = "_4.5sigma_charge_cut"; break;
//     case sigma4_cluster_size_cut: suffix = "_4.5sigma_size_cut"; break;
//     case sigma4_cluster_charge_size_cut: suffix = "_4.5sigma_charge_size_cut"; break;
  }

  // input file
  const bool is_4sigma = 
    cut_id == sigma4_no_cut || cut_id == sigma4_cluster_charge_cut || 
    cut_id == sigma4_cluster_size_cut || cut_id == sigma4_cluster_charge_size_cut;
  const TString inputfilename = is_4sigma ?  Form( "Rootfiles/RawDataClusterTree-%08i-0000-full_4sigma.root", runnumber ) : Form( "Rootfiles/RawDataClusterTree-%08i-0000-full.root", runnumber );
  // const TString inputfilename = is_4sigma ?  Form( "Rootfiles/RawDataClusterTree-%08i-0000-full_4.5sigma.root", runnumber ) : Form( "Rootfiles/RawDataClusterTree-%08i-0000-full.root", runnumber );
  std::cout << "RawDataEfficiency - inputfilename: " << inputfilename << std::endl;
  auto tfile = std::unique_ptr<TFile>( TFile::Open( inputfilename, "READ" ) );
  auto tree = static_cast<TTree*>(tfile->Get("T"));
  
  // output file
  const TString pdffilename = Form( "Figures/RawDataEfficiency_2d_%08i%s.pdf", runnumber, suffix.Data() );
  std::cout << "RawDataEfficiency - pdffilename: " << pdffilename << std::endl;
  PdfDocument pdfdocument( pdffilename );
  
  // loop over detectors
  for( int ilayer = 0; ilayer < 2; ++ilayer )
    for( int itile = 0; itile < 8; ++itile )
  {
    const int detid = itile + 8*ilayer;
    const int detid_ref = (ilayer==0) ? detid+8:detid-8;

    const auto& detname = detector_names[detid];

    const auto var = Form( "first_cluster_strip[%i]", detid_ref );
    
    // reference cut
    TCut ref_cut = Form( "n_detector_clusters[%i]==1", detid_ref );
    if( detid == 5 ) ref_cut = ref_cut && Form( "n_region_clusters[%i]==0", 4*detid_ref+1 );

    TCut first_cluster_strip_cut = Form( "first_cluster_strip[%i]>=100&&first_cluster_strip[%i]<200", detid_ref, detid_ref );
    TCut min_cluster_charge_cut = Form( "min_cluster_charge[%i]>200", detid_ref ); 
    TCut min_cluster_size_cut = Form( "min_cluster_size[%i]>=2", detid_ref );
    
    switch( cut_id )
    {
      default: break;
      
      case cluster_strip_cut:
      ref_cut = ref_cut&&first_cluster_strip_cut;
      break;
      
      case cluster_charge_cut:
      case sigma4_cluster_charge_cut:
      ref_cut = ref_cut&&min_cluster_charge_cut;
      break;
      
      case cluster_size_cut:
      case sigma4_cluster_size_cut:
      ref_cut = ref_cut&&min_cluster_size_cut;
      break;
      
      case cluster_charge_size_cut:
      case sigma4_cluster_charge_size_cut:
      ref_cut = ref_cut&&min_cluster_size_cut&&min_cluster_charge_cut;
      break;
    }
       
    // reference histogram
    const auto hrefname = Form( "h_ref_%s", detname.c_str() );
    auto h_ref = new TH1F( hrefname, "", 128, 0, 256 );
    h_ref->GetXaxis()->SetTitle( "strip" );
    h_ref->GetYaxis()->SetTitle( "counts" );

    Utils::TreeToHisto( tree, h_ref->GetName(), var, ref_cut, false );


    std::cout << "h_ref: " << h_ref << std::endl;
    
    // efficiency cut and found histogram
    const TCut eff_cut = ref_cut && Form( "n_detector_clusters[%i]>0", detid );
    const auto hfoundname = Form( "h_found_%s", detname.c_str() );
    auto h_found = static_cast<TH1*>( h_ref->Clone( hfoundname ) );
    Utils::TreeToHisto( tree, h_found->GetName(), var, eff_cut, false );
    
    // efficiency histogram
    const auto heffname = Form( "h_eff_%s", detname.c_str() );
    auto h_eff = static_cast<TH1*>( h_ref->Clone( heffname ) );
    h_eff->GetYaxis()->SetTitle( "efficiency" );
    Utils::DivideHistograms( h_found, h_ref, h_eff );

    // make plot
    auto cv_name = Form( "cv_%s", detname.c_str() );
    auto cv = new TCanvas( cv_name, "", 1200, 800 );
    cv->Divide( 2, 1 );
    cv->cd(1);
    
    h_ref->SetLineColor(1);
    h_ref->Draw();
    
    h_found->SetLineColor(2);
    h_found->Draw("same");
    Draw::PutText( 0.18, 0.18, detname.c_str());
    
    cv->cd(2);
    h_eff->SetMarkerStyle(20);
    h_eff->SetLineColor(1);
    h_eff->SetMarkerColor(1);
    h_eff->SetMaximum(1);
    h_eff->Draw();
    Draw::PutText( 0.18, 0.18, detname.c_str());
    
    pdfdocument.Add(cv);
    
  }
  
}
