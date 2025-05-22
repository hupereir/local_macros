#include <TFile.h>
#include <TH3.h>

#include <array>
#include <cmath>
  

TString ConvertDistortions()
{

  // this assumes that input bins and output bins are identicals
  
//   const TString inputFile = "/star/u/rcorliss/sphenix/workfest2021/empty_distortion.workfest2021.distortion_map.hist.root";
//   const TString outputFile = "distortion_maps-new/empty_distortion_converted.root";

  const TString inputFile = "/star/u/rcorliss/sphenix/workfest2021/average_minus_static_distortion.workfest2021.distortion_map.hist.root";
  const TString outputFile = "distortion_maps-new/average_minus_static_distortion_converted.root";
  
  auto input = TFile::Open( inputFile );
  
  using hname_pair_t = std::pair<std::string,std::string>;
  const std::vector<hname_pair_t> histogram_names = {
    {"hIntDistortionNegP", "hIntDistortionP_negz"},
    {"hIntDistortionNegR", "hIntDistortionR_negz"},
    {"hIntDistortionNegZ", "hIntDistortionZ_negz"},
    {"hIntDistortionPosP", "hIntDistortionP_posz"},
    {"hIntDistortionPosR", "hIntDistortionR_posz"},
    {"hIntDistortionPosZ", "hIntDistortionZ_posz"}};
    
  auto output = TFile::Open( outputFile, "RECREATE" );
 
  for( const auto& [source,destination]:histogram_names )
  {
    auto hin = dynamic_cast<TH3*>(input->Get(source.c_str()));
    assert( hin );
    hin->SetName(destination.c_str());
    output->cd();
    hin->Write();
  }

  output->Close();
  return outputFile;
}
