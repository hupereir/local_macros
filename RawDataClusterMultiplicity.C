#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH1.h>
#include <TF1.h>

#include <iostream>
#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>

namespace
{
  const std::string logfile = "RawDataClusterMultiplicity.log";
  std::ofstream out( logfile.c_str() );
  bool first = true;
}

//_____________________________________________________________________________
void RawDataClusterMultiplicity(int runNumber = 21162)
{

  set_style( false );

  MicromegasMapping mapping;

  const TString inputFile = Form( "DST/CONDOR_CombinedDataEvaluation/dst_eval_masked-%08i-0000-full.root", runNumber );
  const TString pdfFile = Form( "Figures/RawDataClusterMultiplicity-%08i-0000.pdf", runNumber );
  const TString rootFilename = Form( "Rootfiles/RawDataClusterMultiplicity-%08i-0000.root", runNumber );

  std::cout << "RawDataClusterMultiplicity - inputFile: " << inputFile << std::endl;
  std::cout << "RawDataClusterMultiplicity - pdfFile: " << pdfFile << std::endl;
  std::cout << "RawDataClusterMultiplicity - rootFilename: " << rootFilename << std::endl;
  std::cout << "RawDataClusterMultiplicity - logfile: " << logfile << std::endl;

  PdfDocument pdfDocument( pdfFile );
  RootFile rootfile( rootFilename );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  static constexpr int cluster_max = 50;

  const auto h_phi = new TH1F( "h_phi", "", cluster_max, 0, cluster_max );
  const auto h_z = new TH1F( "h_z", "", cluster_max, 0, cluster_max );
  for( auto&& h:{h_phi,h_z})
  {
    h->StatOverflows(true);
    h->GetXaxis()->SetTitle( "cluster count" );
    h->GetYaxis()->SetTitle( "A.U." );
    h->SetFillStyle(1001);
    h->SetFillColor(kYellow );
    rootfile.Add(h);
  }

  const std::vector<std::string> detectors_phi =
  {
    "SCIP",
    // "NCOP",
    "SEIP",
    "NEIP",
    "SWIP",
    "NWIP"
  };


  const std::vector<std::string> detectors_z =
  {
    "SCIZ",
    // "NCOZ",
    "SEIZ",
    "NEIZ",
    "SWIZ",
    "NWIZ"
  };

  // store mean multiplicity
  std::array<std::string,18> detector_names;
  std::array<double,18> mean_multiplicity_array;

  // create canvas and divide
  const auto cvname = "cv";
  auto cv = new TCanvas( cvname, cvname, 900, 900 );
  cv->Divide(4,4);

  // loop over layers, tiles and region
  for( int ilayer = 0; ilayer < 2; ++ilayer )
  {

    // get actual layer and segmentation
    const int layer = 55+ilayer;
    const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;

    const bool is_phi = (layer == 55);

    for( int itile = 0; itile <8; ++itile )
    {

      const auto idet = itile + 8*ilayer;

      // generate hitset key from layer, tile and segmentation. This allows to retrieve the detector name from the Mapping
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
      detector_names[idet] = name;

      const TString hname( Form( "h_%i", idet ) );
      const TString var = Form( "n_detector_clusters[%i]", idet );
      const auto h = new TH1F( hname, "", cluster_max, 0, cluster_max );
      rootfile.Add(h);
      h->StatOverflows(true);

      TCut full_event_cut( "_events[0]._nrawhits_micromegas>3700" );
      Utils::TreeToHisto( tree, h->GetName(), var, full_event_cut, false );
      h->GetXaxis()->SetTitle( "cluster count" );
      h->GetYaxis()->SetTitle( "A.U." );
      h->SetFillStyle(1001);
      h->SetFillColor(kYellow );

      // sum to average
      if( is_phi )
      {
        if( std::find( detectors_phi.begin(), detectors_phi.end(), name ) != detectors_phi.end() )
        { h_phi->Add( h ); }

      } else {

        if( std::find( detectors_z.begin(), detectors_z.end(), name ) != detectors_z.end() )
        { h_z->Add( h ); }

      }

      cv->cd(idet+1);
      h->Scale( 1./h->GetEntries() );
      h->Draw("hist");

      const double n_cluster_mean = h->GetMean();
      mean_multiplicity_array[idet] = n_cluster_mean;

      gPad->SetLogy( true );

      {
        auto text = new TPaveText(0.2,0.65,0.9,0.84, "NDC" );
        text->SetFillColor(0);
        text->SetFillStyle(0);
        text->SetBorderSize(0);
        text->SetTextAlign(11);
        text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
        text->AddText(Form( "Run %i", runNumber ));
        text->AddText(Form( "%s, #LTclusters#GT=%.1f", name.c_str(), n_cluster_mean ) );
        text->Draw();
      }

    }

  }

  pdfDocument.Add( cv );

  {
    auto cv = new TCanvas( "cvp", "cv", 980, 900 );
    h_phi->Scale( 1./h_phi->GetEntries() );
    h_phi->Draw("hist");

    const double n_cluster_mean = h_phi->GetMean();
    detector_names[16] = "MEAN_P";
    mean_multiplicity_array[16] = n_cluster_mean;

    gPad->SetLogy( true );
    {
      auto text = new TPaveText(0.25,0.71,0.95,0.90, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
      text->AddText(Form( "Run %i", runNumber ));
      text->AddText(Form( "%s, #LTclusters#GT=%.1f", "#Phi strips", n_cluster_mean ) );
      text->Draw();
    }
    pdfDocument.Add(cv);
    // cv->SaveAs(Form( "Figures/RawDataClusterMultiplicity_phi-%08i-0000.pdf", runNumber ) );

  }

  {
    auto cv = new TCanvas( "cvz", "cv", 980, 900 );
    h_z->Scale( 1./h_z->GetEntries() );
    h_z->Draw("hist");

    const double n_cluster_mean = h_z->GetMean();
    detector_names[17] = "MEAN_Z";
    mean_multiplicity_array[17] = n_cluster_mean;
    gPad->SetLogy( true );
    {
      auto text = new TPaveText(0.25,0.71,0.95,0.90, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
      text->AddText(Form( "Run %i", runNumber ));
      text->AddText(Form( "%s, #LTclusters#GT=%.1f", "z strips", n_cluster_mean ) );
      text->Draw();
    }
    pdfDocument.Add(cv);
    // cv->SaveAs(Form( "Figures/RawDataClusterMultiplicity_z-%08i-0000.pdf", runNumber ) );
  }

  // print detectors
  if( first )
  {
    first = false;
    out << "        ";
    for( int i = 0; i < 18; ++i )
    { out << Form( "%8s", detector_names[i].c_str() ); }
    out << std::endl;
  }

  // print occupancies
  out << Form( "%8i", runNumber );
  for( int i = 0; i < 18; ++i )
  { out << Form( "%8.1f", mean_multiplicity_array[i] ); }
  out << std::endl;

  return pdfFile;

}
#include "run_list.h"

//___________________________________
void process_all()
{
  for( const auto& runnumber: RUNS::default_run_list )
  { RawDataClusterMultiplicity( runnumber ); }
}
