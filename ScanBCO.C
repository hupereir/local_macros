#include <RootUtil/FileManager.h>
#include <micromegas/MicromegasRawDataTimingEvaluation.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

void ScanBCO( int runnumber = 45550, int segment = 0)
{

  const TString inputfile = Form( "DST/CONDOR_RawDataEvaluation/MicromegasRawDataTimingEvaluation-%08i-%04i.root", runnumber, segment);
  // const TString inputfile = Form( "DST/CONDOR_RawDataEvaluation/MicromegasRawDataTimingEvaluation-%08i-%04i-test.root", runnumber, segment);
  const TString outputfile = Form( "Rootfiles/ScanBCO-%08i-%04i.root", runnumber, segment);
  std::cout << "ScanBCO - inputfile: " << inputfile << std::endl;
  std::cout << "ScanBCO - outputfile: " << outputfile << std::endl;

  FileManager fileManager( inputfile );
  auto tree = fileManager.GetChain( "T" );
  auto container = new MicromegasRawDataTimingEvaluation::Container;
  tree->SetBranchAddress( "Event", &container );

  // loop over tree entries
  const int entries = tree->GetEntries();
  std::cout << "RawDataClusterTree - entries: " << entries << std::endl;

  bool first = true;
  int n_loop = 0;

  // multiplicator
  const double multiplicator = 4.2629169;
  std::cout << "multiplicator: " << multiplicator << std::endl;

  unsigned int packet_id = 0;

  uint64_t gtm_bco = 0;
  uint64_t gtm_bco_prev = 0;
  uint64_t gtm_bco_first = 0;

  unsigned int fee_bco = 0;
  unsigned int fee_bco_prev = 0;
  unsigned int fee_bco_first = 0;

  uint64_t fee_bco_predicted = 0;

  // output tree
  auto out = new TTree( "out", "out" );
  out->Branch("packet_id", &packet_id);
  out->Branch("gtm_bco", &gtm_bco);
  out->Branch("gtm_bco_prev", &gtm_bco_prev);
  out->Branch("gtm_bco_first", &gtm_bco_first);

  out->Branch("fee_bco", &fee_bco);
  out->Branch("fee_bco_prev", &fee_bco_prev);
  out->Branch("fee_bco_first", &fee_bco_first);
  out->Branch("fee_bco_predicted", &fee_bco_predicted);

  for( int i = 0; i < entries; ++i )
  {

    // some printout
    if( !(i%100) )
    { std::cout << "ScanBCO - entry: " << i << std::endl; }

    tree->GetEntry(i);

    std::cout << "ScanBCO - waveforms: " << container->waveforms.size() << std::endl;
    for( const auto& waveform:container->waveforms )
    {

//       // only look at  one channelm one fee
//       if(
//         !( waveform.fee_id == 9 )&&
//         !( waveform.fee_id == 23 )
//         ) continue;

      // only look at  one channelm one fee
      if(
        !( waveform.fee_id == 9 && waveform.channel == 0 )&&
        !( waveform.fee_id == 23 && waveform.channel == 0 )
        ) continue;

      packet_id = waveform.packet_id;
      fee_bco = waveform.fee_bco;
      gtm_bco = waveform.gtm_bco;
      if( first )
      {
        first = false;

        gtm_bco_prev = gtm_bco;
        gtm_bco_first = gtm_bco;


        fee_bco_first = fee_bco;
        fee_bco_prev = fee_bco;

        std::cout << "scan_bco: gtm_bco_first: " << gtm_bco_first << std::endl;
        std::cout << "scan_bco: fee_bco_first: " << fee_bco_first << std::endl;
        continue;
      }

      fee_bco_predicted = waveform.fee_bco_predicted;
//       std::cout << "ScanBCO -"
//         << " entry: " << i
//         << std::hex
//         << " fee_bco: 0x" << waveform.fee_bco
//         << " fee_bco_predicted: 0x" << fee_bco_predicted
//         << std::dec
//         << std::endl;
//
      out->Fill();
      gtm_bco_prev = gtm_bco;
      fee_bco_prev = fee_bco;
    }

  }

  auto fout = TFile::Open( outputfile, "RECREATE" );
  out->Write();
  fout->Close();

  std::cout << "all done" << std::endl;
}

void run_all()
{
//   for( int i = 2; i < 42; ++i )
//   { ScanBCO(45288,i); }

  for( int runnumber: {45490, 45495, 45550, 45628, 45645 } )
  { ScanBCO( runnumber, 0 ); }

}
