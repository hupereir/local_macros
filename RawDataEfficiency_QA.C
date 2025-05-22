#include <TFile.h>

void RawDataEfficiency_QA()
{

  const TString inputfile = "Rootfiles/QA/MicromegasClusterQA_GL1/MicromegasClusterQA-00053199-0000.root";
  auto f = TFile::Open( inputfile );
  auto href = static_cast<TH1*>( f->Get("h_MicromegasClusterQA_clustercount_ref" ));
  auto hfound = static_cast<TH1*>( f->Get("h_MicromegasClusterQA_clustercount_found" ));

  std::cout << "RawDataEfficiency_QA - inputfile: " << inputfile << std::endl;

  for( int i=0; i<16; ++i )
  {
    const auto detname = href->GetXaxis()->GetBinLabel( i+1 );
    const double ref = href->GetBinContent(i+1);
    const double found = hfound->GetBinContent(i+1);
    const double eff=found/ref;
    const double err = std::sqrt(eff*(1.-eff)/ref);

    std::cout
      << "get_detector_efficiencies -"
      << " detid: " << i
      << " name: " << detname
      << " ref: " << ref
      << " found: " << found
      << " efficiency: " << eff << "+/-" << err
      << std::endl;
  }
}

