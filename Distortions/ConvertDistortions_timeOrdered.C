#include <TFile.h>
#include <TH3.h>

#include <array>
#include <cmath>

namespace
{
  template<class T> inline constexpr T square( const T& x ) { return x*x; }
  template<class T> inline T get_r( const T& x, const T& y ) { return std::sqrt(square(x)+square(y)); }


  // print histogram
  void print_histogram( TH3* h )
  {
    
    std::cout << "ConvertDistortions - name: " << h->GetName() << std::endl;
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
  
TString ConvertDistortions_timeOrdered()
{
  const TString inputFile = TString(getenv("CALIBRATIONROOT")) + "/distortion_maps/TimeOrderedDistortions.root";
  const TString outputFile = "distortion_maps-new/TimeOrderedDistortions_converted.root";
  
  auto input = TFile::Open( inputFile );
  if( !input ) 
  {
    std::cout << "ConvertDistortions_timeOrdered - invalid file: " << inputFile << std::endl;
    return TString();
  }
  
  auto tree = static_cast<TTree*>( input->Get("TimeDists") );
  if( !tree ) 
  {
    std::cout << "ConvertDistortions_timeOrdered - invalid tree" << std::endl;
    return TString();
  }
  
  // create input histograms
  auto hIntDistortionX = new TH3F;
  auto hIntDistortionY = new TH3F;
  auto hIntDistortionZ = new TH3F;
  tree->SetBranchAddress("hIntDistortionX", &hIntDistortionX);
  tree->SetBranchAddress("hIntDistortionY", &hIntDistortionY);
  tree->SetBranchAddress("hIntDistortionZ", &hIntDistortionZ);
 
  // create output TFile, tree and histograms
  auto output = TFile::Open( outputFile, "RECREATE" );
  auto tree_out = new TTree( "TimeDists", "TimeDists" );
  
  std::array<TH3*,2> hDPint;
  std::array<TH3*,2> hDRint;
  std::array<TH3*,2> hDZint;
  
  hDPint[0] = new TH3F( "hIntDistortionP_negz", "hIntDistortionP_negz", 
    82, -0.0785398, 6.36173,
    54, 18.8846, 79.1154,
    82, -106.819, 1.31875 );

  hDRint[0] = new TH3F( "hIntDistortionR_negz", "hIntDistortionR_negz", 
    82, -0.0785398, 6.36173,
    54, 18.8846, 79.1154,
    82, -106.819, 1.31875 );

  hDZint[0] = new TH3F( "hIntDistortionZ_negz", "hIntDistortionZ_negz", 
    82, -0.0785398, 6.36173,
    54, 18.8846, 79.1154,
    82, -106.819, 1.31875 );
  
  hDPint[1] = new TH3F( "hIntDistortionP_posz", "hIntDistortionP_posz", 
    82, -0.0785398, 6.36173,
    54, 18.8846, 79.1154,
    82, -1.31875, 106.819 );

  hDRint[1] = new TH3F( "hIntDistortionR_posz", "hIntDistortionR_posz", 
    82, -0.0785398, 6.36173,
    54, 18.8846, 79.1154,
    82, -1.31875, 106.819 );

  hDZint[1] = new TH3F( "hIntDistortionZ_posz", "hIntDistortionZ_posz", 
    82, -0.0785398, 6.36173,
    54, 18.8846, 79.1154,
    82, -1.31875, 106.819 );
  
  for( int i = 0; i < 2; ++i )
  {

    tree_out->Branch( hDPint[i]->GetName(), "TH3F", &hDPint[i] ); 
    tree_out->Branch( hDRint[i]->GetName(), "TH3F", &hDRint[i] ); 
    tree_out->Branch( hDZint[i]->GetName(), "TH3F", &hDZint[i] ); 
  
  }
  
  // loop over entries
  const auto entries = tree->GetEntries();
  // const int entries = 10;
  std::cout << "ConvertDistortions_timeOrdered - entries: " << entries << std::endl;
  for( int ientry = 0; ientry < entries; ++ientry )
  {
       
    std::cout << "ConvertDistortions_timeOrdered - entry: " << ientry << std::endl;
    tree->GetEntry(ientry);
    
    std::cout << "ConvertDistortions_timeOrdered - x range: " << hIntDistortionX->GetMinimum() << ", " << hIntDistortionX->GetMaximum() << std::endl;
    std::cout << "ConvertDistortions_timeOrdered - y range: " << hIntDistortionY->GetMinimum() << ", " << hIntDistortionY->GetMaximum() << std::endl;
    std::cout << "ConvertDistortions_timeOrdered - z range: " << hIntDistortionZ->GetMinimum() << ", " << hIntDistortionZ->GetMaximum() << std::endl;
    
    // reset all histograms
    for( std::array<TH3*,2> harray: {hDPint, hDRint, hDZint} )
    { 
      harray[0]->Reset();
      harray[1]->Reset();
    }
    
    // fill histograms
    for( int i = 0; i < 2; ++i )
    {
      for( int iz = 0; iz < hDPint[i]->GetNbinsZ(); ++iz )
      {
        
        // get z
        const auto z = hDPint[i]->GetZaxis()->GetBinCenter(iz+1);
        
        // get corrsponding bin in input histogram
        const auto iz_in = hIntDistortionX->GetZaxis()->FindBin( z );
        
        for( int ix = 0; ix < hDPint[i]->GetNbinsX(); ++ix )
          for( int iy = 0; iy < hDPint[i]->GetNbinsY(); ++iy )
        {
          
          // get bin centers
          const auto phi = hDPint[i]->GetXaxis()->GetBinCenter(ix+1);
          const auto r = hDPint[i]->GetYaxis()->GetBinCenter(iy+1);
          
          const auto x = r*std::cos( phi );
          const auto y = r*std::sin( phi );
          
          // get cartesian corrections
          const auto dx = hIntDistortionX->GetBinContent( ix+1, iy+1, iz_in );
          const auto dy = hIntDistortionY->GetBinContent( ix+1, iy+1, iz_in );
          const auto dz = hIntDistortionZ->GetBinContent( ix+1, iy+1, iz_in );
          
          //     const auto dx = hDXint_in->Interpolate( phi, r, std::abs(z) );
          //     const auto dy = hDYint_in->Interpolate( phi, r, std::abs(z) );
          //     const auto dz = hDZint_in->Interpolate( phi, r, std::abs(z) );
          
          
          // calculate polar distortions
          const auto dr = get_r( x+dx, y+dy ) - get_r( x, y );
          auto dphi = std::atan2( y+dy, x+dx ) - std::atan2( y, x );
          if( dphi < -M_PI ) dphi += 2*M_PI;
          if( dphi >= M_PI ) dphi -= 2*M_PI;
          
          // fill histograms
          hDPint[i]->SetBinContent( ix+1, iy+1, iz+1, r*dphi );
          hDRint[i]->SetBinContent( ix+1, iy+1, iz+1, dr );
          hDZint[i]->SetBinContent( ix+1, iy+1, iz+1, dz );    
        }
        
      }
      
    }
  
    // save in tree
    tree_out->Fill();
  }
    
  // write to output
  output->cd();
  tree_out->Write();
  
  return outputFile;

}
