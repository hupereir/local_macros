#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH1.h>
#include <TF1.h>

#include <memory>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>

namespace
{
  const std::string logfile = "RawDataClusterCharge.log";
  std::ofstream out( logfile.c_str() );
  bool first = true;
}

//_____________________________________________________________________________
void RawDataClusterCharge(int runNumber = 21162)
{

  set_style( false );

  const TString inputFile = Form( "DST/CONDOR_CombinedDataEvaluation/dst_eval_masked-%08i-0000-full.root", runNumber );
  const TString pdfFile = Form( "Figures/RawDataClusterCharge-%08i-0000.pdf", runNumber );
  const TString rootFilename = Form( "Rootfiles/RawDataClusterCharge-%08i-0000.root", runNumber );


  std::cout << "RawDataClusterCharge - inputFile: " << inputFile << std::endl;
  std::cout << "RawDataClusterCharge - pdfFile: " << pdfFile << std::endl;
  std::cout << "RawDataClusterCharge - rootFilename: " << rootFilename << std::endl;
  std::cout << "RawDataClusterCharge - logfile: " << logfile << std::endl;

  PdfDocument pdfDocument( pdfFile );
  RootFile rootfile( rootFilename );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  MicromegasMapping mapping;

  // store mean cluster charge
  std::array<std::string,18> detector_names;
  std::array<double,18> mean_charge_array;

  // create canvas and divide
  const auto cvname = "cv";
  auto cv = new TCanvas( cvname, cvname, 900, 900 );
  cv->Divide(4,4);

  static constexpr double charge_max = 2e3;
  static constexpr int ncharge_bins = 100;
  static constexpr bool use_log = false;

  const auto h_phi = new TH1F( "h_phi", "", ncharge_bins, 0, charge_max );
  const auto h_z = new TH1F( "h_z", "", ncharge_bins, 0, charge_max );
  for( auto&& h:{h_phi,h_z})
  {

    h->StatOverflows(true);

    h->GetXaxis()->SetTitle( "cluster charge (ADC)" );
    h->GetYaxis()->SetTitle( "A.U." );
    h->SetFillStyle(1001);
    h->SetFillColor(kYellow );
    rootfile.Add(h);
  }

  const std::vector<std::string> detectors_phi =
  {
    "SCIP",
    "NCOP",
    "SEIP",
    "NEIP",
    "SWIP",
    "NWIP"
  };


  const std::vector<std::string> detectors_z =
  {
    "SCIZ",
    "NCOZ",
    "SEIZ",
    "NEIZ",
    "SWIZ",
    "NWIZ"
  };

  TGaxis::SetMaxDigits(4);

  // loop over layers, tiles and region
  for( int ilayer = 0; ilayer < 2; ++ilayer )
  {

    // get actual layer and segmentation
    const int layer = 55+ilayer;
    const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;

    const bool is_phi = (layer == 55);

    for( int itile = 0; itile <8; ++itile )
    {

      // generate hitset key from layer, tile and segmentation. This allows to retrieve the detector name from the Mapping
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );

      std::cout << "RawDataClusterCharge - name: " << name << std::endl;

      const auto idet = itile + 8*ilayer;
      const int bin = idet;
      detector_names[idet] = name;

      const auto hname = Form( "h2d_%i", idet );
      auto h = new TH1F( hname, "", ncharge_bins, 0, charge_max );
      h->StatOverflows(true);
      h->GetXaxis()->SetTitle( "cluster charge (ADC)" );
      h->GetYaxis()->SetTitle( "A.U." );
      h->SetFillStyle(1001);
      h->SetFillColor(kYellow );
      rootfile.Add(h);

      const TString var = "clusters[].charge";
      TCut det_cut = Form( "(clusters[].tile+8*(clusters[].layer-55))==%i", idet );
      TCut full_event_cut( "_events[0]._nrawhits_micromegas>3700" );
      Utils::TreeToHisto( tree, h->GetName(), var, full_event_cut&&det_cut, false );

      // sum to average
      if( is_phi )
      {
        if( std::find( detectors_phi.begin(), detectors_phi.end(), name ) != detectors_phi.end() )
        { h_phi->Add( h ); }

      } else {

        if( std::find( detectors_z.begin(), detectors_z.end(), name ) != detectors_z.end() )
        { h_z->Add( h ); }

      }

      cv->cd(bin+1);
      h->Scale( 1./h->GetEntries() );
      h->Draw("hist");

      const double charge_mean = h->GetMean();
      mean_charge_array[idet] = charge_mean;

      if( use_log ) gPad->SetLogy( true );
      gPad->SetRightMargin( 0.10 );

      {
        auto text = new TPaveText(0.37,0.76,0.94,0.9, "NDC" );
        text->SetFillColor(0);
        text->SetFillStyle(0);
        text->SetBorderSize(0);
        text->SetTextAlign(11);
        text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
        text->AddText(Form( "Run %i", runNumber ));
        text->AddText(Form( "%s, #LTcluster charge#GT=%.1f", name.c_str(), charge_mean ) );
        text->Draw();
      }

    }

  }

  pdfDocument.Add( cv );

  {
    auto cv = new TCanvas( "cvp", "cv", 980, 900 );
    h_phi->Scale( 1./h_phi->GetEntries() );
    h_phi->Draw("hist");

    const double charge_mean = h_phi->GetMean();
    detector_names[16] = "MEAN_P";
    mean_charge_array[16] = charge_mean;
    if( use_log ) gPad->SetLogy( true );
    gPad->SetRightMargin( 0.10 );
    {
      auto text = new TPaveText(0.37,0.76,0.94,0.9, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
      text->AddText(Form( "Run %i", runNumber ));
      text->AddText( Form("#Phi strips, #LTcluster charge#GT=%.1f", charge_mean ) );
      text->Draw();
    }
    pdfDocument.Add(cv);
    cv->SaveAs(Form( "Figures/RawDataClusterCharge_phi-%08i-0000.pdf", runNumber ) );

  }

  {
    auto cv = new TCanvas( "cvz", "cv", 980, 900 );
    h_z->Scale( 1./h_z->GetEntries() );
    h_z->Draw("hist");

    const double charge_mean = h_z->GetMean();
    detector_names[17] = "MEAN_Z";
    mean_charge_array[17] = charge_mean;

    if( use_log ) gPad->SetLogy( true );
    gPad->SetRightMargin( 0.10 );
    {
      auto text = new TPaveText(0.37,0.76,0.94,0.9, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
      text->AddText(Form( "Run %i", runNumber ));
      text->AddText(Form( "z strips, #LTcluster charge#GT=%.1f", charge_mean ) );
      text->Draw();
    }
    pdfDocument.Add(cv);
    cv->SaveAs(Form( "Figures/RawDataClusterCharge_z-%08i-0000.pdf", runNumber ) );
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
  { out << Form( "%8.1f", mean_charge_array[i] ); }
  out << std::endl;

  return pdfFile;

}

#include "run_list.h"

//___________________________________
void process_all()
{
  for( const auto& runnumber: RUNS::default_run_list )
  { RawDataClusterCharge( runnumber ); }
}
