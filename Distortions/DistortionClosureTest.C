
#include <TFile.h>
#include <TH3.h>

#include <array>
#include <cmath>

namespace
{
  // print histogram
  void print_histogram( TH3* h )
  {
    
    std::cout << "InvertDistortions - name: " << h->GetName() << std::endl;
    for( const auto& axis:{h->GetXaxis(), h->GetYaxis(), h->GetZaxis() } )
    {
      std::cout 
        << "  " << axis->GetName() 
        << " bins: " << axis->GetNbins() 
        << " min: " << axis->GetXmin() 
        << " max: " << axis->GetXmax() 
        << std::endl;
    }
    std::cout << std::endl;
  }
  
}

TString DistortionClosureTest(
  const TString inputFile = "/star/u/rcorliss/sphenix/trackingStudySampleNov2021/static_only.distortion_map.hist.root",
//   const TString inputFile_inverted = "/star/u/rcorliss/sphenix/trackingStudySampleNov2021/static_only.distortion_map.hist.root",
//   const TString outputFile = "distortion_maps-new/static_only_closure.root"
  const TString inputFile_inverted = "distortion_maps-new/static_only_inverted_4-new.root",
  const TString outputFile = "distortion_maps-new/static_only_closure.root"
  
//   const TString inputFile = "distortion_maps-new/average_minus_static_distortion_converted.root",
//   const TString inputFile_inverted = "distortion_maps-new/average_minus_static_distortion_inverted_4-new.root",
//   const TString outputFile = "distortion_maps-new/average_minus_static_distortion_closure_4-new.root"
  )
{
  
  std::cout << "InvertDistortionMap - inputFile: " << inputFile << std::endl;
  std::cout << "InvertDistortionMap - inputFile_inverted: " << inputFile_inverted << std::endl;
  std::cout << "InvertDistortionMap - outputFile: " << outputFile << std::endl;
  
  auto input = TFile::Open( inputFile );
  auto input_inverted = TFile::Open( inputFile_inverted );
  auto output = TFile::Open( outputFile, "RECREATE" );
        
  auto hdp_closure_1d = new TH1F( "dp_closure_1d", "dp_closure_1d", 100, -0.01, 0.01 );
  auto hdr_closure_1d = new TH1F( "dr_closure_1d", "dr_closure_1d", 100, -0.01, 0.01 );
  auto hdz_closure_1d = new TH1F( "dz_closure_1d", "dz_closure_1d", 100, -0.01, 0.01 );
  
  std::array<TString,2> suffix = { "negz", "posz" };
  for( int iside = 0; iside < 2; ++iside )
  {
    
    // load relevan histograms from inputFile
    auto hdp = static_cast<TH3*>( input->Get( Form( "hIntDistortionP_%s", suffix[iside].Data() ) ) );
    auto hdr = static_cast<TH3*>( input->Get( Form( "hIntDistortionR_%s", suffix[iside].Data() ) ) );
    auto hdz = static_cast<TH3*>( input->Get( Form( "hIntDistortionZ_%s", suffix[iside].Data() ) ) );
    print_histogram( hdp );
    
    // load relevan histograms from inputFile
    auto hdp_inverted = static_cast<TH3*>( input->Get( Form( "hIntDistortionP_%s", suffix[iside].Data() ) ) );
    auto hdr_inverted = static_cast<TH3*>( input_inverted->Get( Form( "hIntDistortionR_%s", suffix[iside].Data() ) ) );
    auto hdz_inverted = static_cast<TH3*>( input_inverted->Get( Form( "hIntDistortionZ_%s", suffix[iside].Data() ) ) );

    // create closure histograms.
    auto hdp_closure = static_cast<TH3*>( hdp->Clone( Form( "hIntDistortionP_%s_closure", suffix[iside].Data() ) ) );
    auto hdr_closure = static_cast<TH3*>( hdp->Clone( Form( "hIntDistortionR_%s_closure", suffix[iside].Data() ) ) );
    auto hdz_closure = static_cast<TH3*>( hdp->Clone( Form( "hIntDistortionZ_%s_closure", suffix[iside].Data() ) ) );
    
    // reset
    for( const auto& h:{hdp_closure, hdr_closure, hdz_closure} ) { h->Reset(); }

    // get then number and limits of active bins
    const int npbins = hdp->GetNbinsX()-2;
    const int nrbins = hdp->GetNbinsY()-2;
    const int nzbins = hdp->GetNbinsZ()-2;
    
    // scan inputs histograms
    for( int ip = 0; ip < npbins; ++ip )
    {
      for( int ir = 0; ir < nrbins; ++ir )
      {
        for( int iz = 0; iz < nzbins; ++iz )
        {
          
          // get initial position
          const double p = hdp->GetXaxis()->GetBinCenter( ip+2 );
          const double r = hdp->GetYaxis()->GetBinCenter( ir+2 );
          const double z = hdp->GetZaxis()->GetBinCenter( iz+2 );
          
          //  get distortions from source histograms
          const double drp =  hdp->Interpolate(p,r,z);
          const double dp = drp/r;
          const double dr = hdr->Interpolate(p,r,z);
          const double dz = hdz->Interpolate(p,r,z);
          
          // get new position
          double p_new = p+dp;
          if( p_new < 0 ) p_new += 2.*M_PI;
          if( p_new >= 2*M_PI ) p_new -= 2.*M_PI;
          
          const double r_new = r+dr;
          const double z_new = z+dz;

          // get inverted distortions
          const double drp_inverted =  hdp_inverted->Interpolate(p_new,r_new,z_new);
          const double dp_inverted = drp_inverted/r_new;
          const double dr_inverted = hdr_inverted->Interpolate(p_new,r_new,z_new);
          const double dz_inverted = hdz_inverted->Interpolate(p_new,r_new,z_new);
                 
          // update positions
          double p_closure = p_new - dp_inverted;
          if( p_closure < 0 ) p_closure += 2.*M_PI;
          if( p_closure >= 2*M_PI ) p_closure -= 2.*M_PI;
          const double r_closure = r_new - dr_inverted;
          const double z_closure = z_new - dz_inverted;
          
          // calculate closure
          double dp_closure = p_closure - p;
          if( dp_closure < -M_PI ) dp_closure += 2.*M_PI;
          if( dp_closure >= M_PI ) dp_closure -= 2.*M_PI;
          const double dr_closure = r_closure - r;
          const double dz_closure = z_closure -z;
          
          // fill converted histogram
          hdp_closure->Fill( p, r, z, dp_closure*r );
          hdr_closure->Fill( p, r, z, dr_closure );
          hdz_closure->Fill( p, r, z, dz_closure );

          // also fill 1d histograms
          hdp_closure_1d->Fill( dp_closure*r );
          hdr_closure_1d->Fill( dr_closure );
          hdz_closure_1d->Fill( dz_closure );
                    
        }
      }
    }

    // write to ouput TFile
    output->cd();
    for( const auto& h:{hdp_closure, hdr_closure, hdz_closure} ) 
    { h->Write(); }
    
  }

  // write 1D histograms to ouput TFile, using new name
  output->cd();
  for( const auto& h:{hdp_closure_1d, hdr_closure_1d, hdz_closure_1d} ) 
  { h->Write(); }
  
  output->Close();
  return outputFile;
}
