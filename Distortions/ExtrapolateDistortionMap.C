#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH3.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

// define reference phi r and z for the extrapolation
static constexpr double isec_ref = 3;
static constexpr double phi_ref = isec_ref*M_PI/6 + M_PI/12;

// radius of the micromegas layer
static constexpr double r_ref = 82;

// z extrapolation between micromegas
// static constexpr double zextrap_min = 51.25;
// static constexpr double zextrap_max = 53.75;
static constexpr double zextrap_min = 48;
static constexpr double zextrap_max = 58;

// Micromegas acceptance in incomplete sectors
static constexpr double zref = 33.25;
static constexpr double length = 50 - 5;
static constexpr double zref_min = zref - length/2;
static constexpr double zref_max = zref + length/2;


namespace 
{
  template< class T>
    inline constexpr T square( const T& x ) { return x*x; }
}

//_______________________________________________
/* z extrapolation
 * interpolate between micromegas in the fully equiped sector
 */
void ExtrapolateZ( TH3* hin )
{
  if( !hin ) return;

  Utils::PrintAxis( hin );  

  // get reference phi bin
  const int phibin_ref = hin->GetXaxis()->FindBin( phi_ref );

  // loop over radial bins
  for( int ir = 0; ir < hin->GetYaxis()->GetNbins(); ++ir )
  {
    
    // get current radius
    const auto r = hin->GetYaxis()->GetBinCenter( ir+1 );
    
    // get z integration window for reference
    const auto zextrap_min_loc = zextrap_min * r/r_ref;
    const auto zextrap_max_loc = zextrap_max * r/r_ref;
  
    // get corresponding bins
    const int zbin_min[2] = { hin->GetZaxis()->FindBin( -zextrap_max_loc ), hin->GetZaxis()->FindBin( zextrap_min_loc ) };
    const int zbin_max[2] = { hin->GetZaxis()->FindBin( -zextrap_min_loc ), hin->GetZaxis()->FindBin( zextrap_max_loc ) };
    
    for( int isign = 0; isign < 2; ++isign )
    {
      // adjust z positions
      const auto z_min = hin->GetZaxis()->GetBinCenter( zbin_min[isign] );
      const auto z_max = hin->GetZaxis()->GetBinCenter( zbin_max[isign] );
    
      // get reference
      const auto content_min = hin->GetBinContent( phibin_ref, ir+1, zbin_min[isign] );
      const auto content_max = hin->GetBinContent( phibin_ref, ir+1, zbin_max[isign] );
      const auto error_min = hin->GetBinError( phibin_ref, ir+1, zbin_min[isign] );
      const auto error_max = hin->GetBinError( phibin_ref, ir+1, zbin_max[isign] );

      // loop over z bins
      for( int iz = zbin_min[isign]+1; iz < zbin_max[isign]; ++iz )
      {
      
        const auto z = hin->GetZaxis()->GetBinCenter( iz );
        
        // interpolate
        const auto alpha_min = (z_max-z)/(z_max-z_min);
        const auto alpha_max = (z-z_min)/(z_max-z_min);
        
        const auto content = alpha_min*content_min + alpha_max*content_max;
        const auto error = std::sqrt(square( alpha_min * error_min ) + square( alpha_max*error_max));
        
        hin->SetBinContent( phibin_ref, ir+1, iz, content );
        hin->SetBinError( phibin_ref, ir+1, iz, error );
      }
    }
  }
  return;
}
    
//_______________________________________________
/* first phi extrapolation
 * copy the full z dependence of reference sector to all other sectors, separately for positive and negative z, 
 * normalized by the measurement from provided micromegas, at the appropriate z
 */
void ExtrapolatePhi1( TH3* hin )
{
  if( !hin ) return;
  Utils::PrintAxis( hin );  

  // get reference phi bin
  const int phibin_ref = hin->GetXaxis()->FindBin( phi_ref );
  
  // loop over sectors
  for( int isec = 0; isec < 12; ++isec )
  {
        
    // skip reference sector
    if( isec == isec_ref ) continue;
    
    // get relevant phi and corresponding bin
    const double phi = isec*M_PI/6 + M_PI/12;
    const int phibin = hin->GetXaxis()->FindBin( phi );
    
    // loop over radial bins
    for( int ir = 0; ir < hin->GetYaxis()->GetNbins(); ++ir )
    {
      
      // get current radius
      const auto r = hin->GetYaxis()->GetBinCenter( ir+1 );
      
      // get z integration window for reference
      const auto zref_min_loc = zref_min * r/r_ref;
      const auto zref_max_loc = zref_max * r/r_ref;
      
      // get corresponding bins
      const int zbin_ref_neg[2] = { hin->GetZaxis()->FindBin( -zref_max_loc ), hin->GetZaxis()->FindBin( -zref_min_loc ) };      
      const int zbin_ref_pos[2] = { hin->GetZaxis()->FindBin( zref_min_loc ), hin->GetZaxis()->FindBin( zref_max_loc ) };
          
      // loop over z bins
      for( int iz = 0; iz < hin->GetZaxis()->GetNbins(); ++iz )
      {
        const auto content_ref = hin->GetBinContent( phibin_ref, ir+1, iz+1 );
        const auto error_ref = hin->GetBinError( phibin_ref, ir+1, iz+1 );

        #if true 
        // calculate scale factor
        const auto z = hin->GetZaxis()->GetBinCenter( iz+1 );
        const auto norm_ref = hin->Integral( phibin_ref, phibin_ref, ir+1, ir+1, (z>0) ? zbin_ref_pos[0]:zbin_ref_neg[0], (z>0) ? zbin_ref_pos[1]:zbin_ref_neg[1] );  
        const auto norm_loc = hin->Integral( phibin, phibin, ir+1, ir+1, (z>0) ? zbin_ref_pos[0]:zbin_ref_neg[0], (z>0) ? zbin_ref_pos[1]:zbin_ref_neg[1] );  
        const auto scale = (norm_ref == 0) ? 1:norm_loc/norm_ref;
        #else
        const auto scale = 1;
        #endif
        
        // assign to output histogram
        hin->SetBinContent( phibin, ir+1, iz+1, content_ref*scale );
        hin->SetBinError( phibin, ir+1, iz+1, error_ref*scale );        
      }
    }
  }
    
  return; 
}

//_______________________________________________
/* 
 * second phi extrapolation 
 * for each r, z and phi bin, linearly extrapolate between neighbor phi sector measurements 
 */
void ExtrapolatePhi2( TH3* hin )
{
  if( !hin ) return;
  
  for( int iphi = 0; iphi < hin->GetXaxis()->GetNbins(); ++iphi )
  {

    // find nearest sector phi bins
    const auto phi = hin->GetXaxis()->GetBinCenter( iphi+1 );
    const int isec = std::floor( (phi - M_PI/12)/(M_PI/6) );
    double phi_min =  isec*M_PI/6 + M_PI/12;
    double phi_max =  phi_min + M_PI/6;
    
    if( phi_min < 0 ) phi_min += 2*M_PI;
    if( phi_max >= 2*M_PI ) phi_max -= 2*M_PI;
    
    const auto phibin_min = hin->GetXaxis()->FindBin( phi_min );
    if( phibin_min == iphi+1 ) continue;

    const auto phibin_max = hin->GetXaxis()->FindBin( phi_max );    
    if( phibin_max == iphi+1 ) continue;
    
    // loop over radial bins
    for( int ir = 0; ir < hin->GetYaxis()->GetNbins(); ++ir )
    {
      
      // loop over z bins
      for( int iz = 0; iz < hin->GetZaxis()->GetNbins(); ++iz )
      {
        const auto content_min = hin->GetBinContent( phibin_min, ir+1, iz+1 );
        const auto content_max = hin->GetBinContent( phibin_max, ir+1, iz+1 );
        const auto error_min = hin->GetBinError( phibin_min, ir+1, iz+1 );
        const auto error_max = hin->GetBinError( phibin_max, ir+1, iz+1 );
        
        // perform linear extrapolation
        const auto alpha_min = (phi_max-phi)/(phi_max-phi_min);
        const auto alpha_max = (phi-phi_min)/(phi_max-phi_min);
        
        const auto content = alpha_min*content_min + alpha_max*content_max;
        const auto error = std::sqrt(square( alpha_min * error_min ) + square( alpha_max*error_max));
        
        hin->SetBinContent( iphi+1, ir+1, iz+1, content );
        hin->SetBinError( iphi+1, ir+1, iz+1, error );
      }
        
    }
  }
    
  return;
}

//_______________________________________________
TString ExtrapolateDistortionMap()
{
  set_style( false );

  // input grid
  // const TString tag = "_full_realistic_micromegas_mm-new";
  const TString tag = "_full_realistic_micromegas_mm-coarse-new2";
  const auto inputfilename = Form( "Rootfiles/Distortions%s.root", tag.Data() );
  const auto outputfilename =  Form( "Rootfiles/Distortions%s_extrapolated.root", tag.Data() );
  
  std::cout << "ExtrapolateDistortionMap - inputfilename: " << inputfilename << std::endl;
  std::cout << "ExtrapolateDistortionMap - outputfilename: " << outputfilename << std::endl;
  
  auto f = TFile::Open( inputfilename );
  if( !f ) return TString();

  // load input histogram  
  #if true
  auto hentries = dynamic_cast<TH3*>(f->Get("hentries_rec"));
  auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hDistortionP_rec"));
  auto hDistortionR_rec = dynamic_cast<TH3*>(f->Get("hDistortionR_rec"));
  auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hDistortionZ_rec"));
  #else
  TH3* hentries = nullptr;
  auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionP"));
  auto hDistortionR_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionR"));
  auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionZ"));
  #endif
  
  // output grid
  RootFile fout( outputfilename );

  for( const auto& h:{hentries, hDistortionP_rec, hDistortionZ_rec, hDistortionR_rec } )
  {
    if( !h ) continue;
    ExtrapolateZ( h );
    ExtrapolatePhi1( h );
    ExtrapolatePhi2( h );
    fout.Add( h );    
  }
    
  return outputfilename;
}
