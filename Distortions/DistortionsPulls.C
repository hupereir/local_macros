#include <RootUtil/FileManager.h>
#include <RootUtil/RootFile.h>
#include <g4eval/TrackingEvaluator_hp.h>

#include <TTree.h>

R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

#include "LayerDefines.h"

static constexpr bool use_truth_information = false;
static constexpr bool use_micromegas = true;

// histogram limits
static constexpr float max_alpha = 0.5;
static constexpr float max_beta = 1.0;

// bins
int m_phibins = 0;
int m_rbins = 0;
int m_zbins = 0;
int m_totalbins = 0;

//_________________________________________________________________________
template< class T>
constexpr T square( const T& x ) { return x*x; }

//_________________________________________________________________________
template< class T>
T delta_phi( const T& phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//__________________________________________________________________
bool accept_track( const TrackingEvaluator_hp::TrackStruct& track )
{
  
  // momentum cut
  if( track._truth_pt < 0.5 ) return false;
  if( track._pt < 0.5 ) return false;
  
  // hit pattern cuts
  if( track._nclusters_mvtx < 3 ) return false;
  if( track._nclusters_intt < 2 ) return false;
  if( use_micromegas && track._nclusters_micromegas < 2 ) return false;
  
  // also cut on the ndf
  if( use_micromegas && track._ndf < 9 ) return false;

  return true;
}

//__________________________________________________________________
bool accept_cluster( const TrackingEvaluator_hp::ClusterStruct& cluster )
{

  // skip clusters not in TPC
  if( cluster._layer < firstLayer_tpc || cluster._layer >= firstLayer_tpc + nLayers_tpc )
  { return false; }
  
  if( true )
  {
    /*
    remove clusters with too small errors since they are likely pathological
    and have a large contribution to the chisquare
    */
    if( cluster._trk_r*cluster._phi_error < 0.015 ) return false;
    if( cluster._z_error < 0.05 ) return false;
  }

  // check against window
  if( use_truth_information )
  {
    if( std::abs( cluster._truth_alpha ) > max_alpha ) return false;
    if( std::abs( cluster._truth_beta ) > max_beta ) return false;
  } else {
    if( std::abs( cluster._trk_alpha ) > max_alpha ) return false;
    if( std::abs( cluster._trk_beta ) > max_beta ) return false;
  }

  return true;    
}

//_________________________________________________________________________
int get_cell( int iphi, int ir, int iz )
{
  if( iphi < 0 || iphi >= m_phibins ) return -1;
  if( ir < 0 || ir >= m_rbins ) return -1;
  if( iz < 0 || iz >= m_zbins ) return -1;
  return iz + m_zbins*( ir + m_rbins*iphi );
}

//_________________________________________________________________________
int get_cell( TH3* h, const TrackingEvaluator_hp::ClusterStruct& cluster )
{
  
  // phi range
  constexpr float m_phimin = 0;
  constexpr float m_phimax = 2.*M_PI;

  // phi
  auto phi = cluster._phi;
  while( phi < m_phimin ) phi += 2.*M_PI;
  while( phi >= m_phimax ) phi -= 2.*M_PI;
  
  auto iphi = h->GetXaxis()->FindBin( phi );
  if( iphi == 0 || iphi > h->GetNbinsX() ) return -1;
  
  // radius
  auto ir = h->GetYaxis()->FindBin( cluster._r );
  if( ir == 0 || ir > h->GetNbinsY() ) return -1;
  
  // z
  auto iz = h->GetZaxis()->FindBin( cluster._z );
  if( iz == 0 || iz > h->GetNbinsZ() ) return -1;
  
  return get_cell( iphi-1, ir-1, iz-1 );
}

//_________________________________________________________________________
double get_correction( TH3* h, const TrackingEvaluator_hp::ClusterStruct& cluster )
{
  
  // phi range
  constexpr float m_phimin = 0;
  constexpr float m_phimax = 2.*M_PI;

  // phi
  auto phi = cluster._phi;
  while( phi < m_phimin ) phi += 2.*M_PI;
  while( phi >= m_phimax ) phi -= 2.*M_PI;
  
  auto iphi = h->GetXaxis()->FindBin( phi );
  if( iphi == 0 || iphi > h->GetNbinsX() ) return 0;
  
  // radius
  auto ir = h->GetYaxis()->FindBin( cluster._r );
  if( ir == 0 || ir > h->GetNbinsY() ) return 0;
  
  // z
  auto iz = h->GetZaxis()->FindBin( cluster._z );
  if( iz == 0 || iz > h->GetNbinsZ() ) return 0;
  
  return h->GetBinContent( iphi, ir, iz );
}

//_________________________________________________________________________
TString DistortionsPulls()
{

  // input files
  const TString tag = "_realistic_micromegas";
  const TString inputFile = "DST/CONDOR_realistic_micromegas/dst_reco_truth_notpc_distortions/dst_reco_realistic_micromegas_*.root";
  
  // root file
  TString correctionFilename;
  if( use_truth_information ) correctionFilename = Form( "Rootfiles/Distortions%s_truth.root", tag.Data() );
  else if( use_micromegas ) correctionFilename = Form( "Rootfiles/Distortions%s_mm.root", tag.Data() );
  else correctionFilename = Form( "Rootfiles/Distortions%s_all.root", tag.Data() );
  auto f = TFile::Open( correctionFilename );
  if( !f ) return TString();

  // output filename
  TString outputFilename;
  if( use_truth_information ) outputFilename = Form( "Rootfiles/DistortionsPulls%s_truth.root", tag.Data() );
  else if( use_micromegas ) outputFilename = Form( "Rootfiles/DistortionsPulls%s_mm.root", tag.Data() );
  else outputFilename = Form( "Rootfiles/DistortionsPulls%s_all.root", tag.Data() );
  
  // load histograms and get dimensions
  auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hDistortionP_rec"));
  auto hDistortionR_rec = dynamic_cast<TH3*>(f->Get("hDistortionR_rec"));
  auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hDistortionZ_rec"));
  if( !( hDistortionP_rec && hDistortionR_rec && hDistortionZ_rec ) ) return TString();
  
  // get axis dimentions
  m_phibins = hDistortionP_rec->GetXaxis()->GetNbins();
  m_rbins = hDistortionP_rec->GetYaxis()->GetNbins();
  m_zbins = hDistortionP_rec->GetZaxis()->GetNbins();
  m_totalbins = m_zbins*m_phibins*m_rbins;
  std::cout << "DistortionsPulls - bins: (" << m_phibins << ", " << m_rbins << ", " << m_zbins << ")" << std::endl;
  
  // global histograms, for quick check
  auto pulls_drphi_all = new TH1F( "pulls_drphi_all", "pulls_drphi_all", 100, -10, 10 );
  pulls_drphi_all->GetXaxis()->SetTitle( "r.#Delta#phi_{cluster-track} pulls" );
    
  auto pulls_dz_all = new TH1F( "pulls_dz_all", "pulls_dz_all", 100, -10, 10 );
  pulls_dz_all->GetXaxis()->SetTitle( "#Deltaz_{cluster-track} pulls" );

  // histogram arrays
  using TH1_ptr_t = std::unique_ptr<TH1>;
  std::vector<TH1_ptr_t> pulls_drphi( m_totalbins );
  std::vector<TH1_ptr_t> pulls_dz( m_totalbins );

  for( int iphi = 0; iphi < m_phibins; ++iphi )
    for( int ir = 0; ir < m_rbins; ++ir )
    for( int iz = 0; iz < m_zbins; ++iz )
  {
    const auto icell = get_cell( iphi, ir, iz );
    
    {
      // rphi pulls
      const auto hname = Form( "pulls_drphi_p%i_r%i_z%i", iphi, ir, iz );
      pulls_drphi[icell].reset( new TH1F( hname, hname, 100, -10, 10 ) );
      pulls_drphi[icell] ->GetXaxis()->SetTitle( "r.#Delta#phi_{cluster-track} pulls" );
    }
    
    {
      // z pulls
      const auto hname = Form( "pulls_dz_p%i_r%i_z%i", iphi, ir, iz );
      pulls_dz[icell].reset( new TH1F( hname, hname, 100, -10, 10 ) );
      pulls_dz[icell]->GetXaxis()->SetTitle( "#Deltaz_{cluster-track} pulls" );
    }

  }

  // get tree
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return TString();

  auto container = new TrackingEvaluator_hp::Container;
  tree->SetBranchAddress( "DST#EVAL#TrackingEvaluator_hp::Container", &container );

  // loop over entries
  const int entries = tree->GetEntries();
  // const int entries = 50000;
  std::cout << "DistortionsPulls - entries: " << entries << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    tree->GetEntry(i);
    if( !(i%100) )
      std::cout << "DistortionsPulls - entry: " << i << std::endl;

    // loop over tracks
    for( const auto& track:container->tracks() )
    {

      // check
      if( !accept_track( track ) ) continue;
      
      // loop over clusters
      for( const auto& cluster:track._clusters )
      {
        
        if( !accept_cluster( cluster ) ) continue;

        // store errors
        const auto erp = use_truth_information ?
          cluster._trk_r*cluster._phi_error:
          cluster._trk_r*std::sqrt( square( cluster._phi_error ) + square( cluster._trk_phi_error ) );
        
        const auto ez = use_truth_information ? 
          cluster._z_error:
          std::sqrt( square( cluster._z_error ) + square( cluster._trk_z_error ) );

        // sanity check
        if( std::isnan( erp ) ) continue;
        if( std::isnan( ez ) ) continue;

        const auto drp = use_truth_information ?
          cluster._truth_r*( delta_phi( cluster._phi - cluster._truth_phi ) ):
          cluster._trk_r*( delta_phi( cluster._phi - cluster._trk_phi ) );
        if( std::isnan(drp) ) continue;

        const auto talpha = use_truth_information ?
          -std::tan( cluster._truth_alpha ):
          -std::tan( cluster._trk_alpha );
        if( std::isnan(talpha) ) continue;

        const auto dz = use_truth_information ? (cluster._z - cluster._truth_z):(cluster._z - cluster._trk_z );
        if( std::isnan(dz) ) continue;

        const auto tbeta = use_truth_information ?
          -std::tan( cluster._truth_beta ):
          -std::tan( cluster._trk_beta );
        if( std::isnan(tbeta) ) continue;

        // find corrections
        const auto drp_corr = get_correction( hDistortionP_rec, cluster );
        const auto dr_corr = get_correction( hDistortionR_rec, cluster );
        const auto dz_corr = get_correction( hDistortionZ_rec, cluster );
        
        // get relevant cell
        const auto icell = get_cell( hDistortionP_rec, cluster );
        if( icell < 0 ) continue;
        
        // fill normalized residuals
        {
          const auto pull = (drp - (drp_corr + talpha*dr_corr))/erp;
          pulls_drphi_all->Fill(pull);
          pulls_drphi[icell]->Fill(pull);
        }
        
        {
          const auto pull = (dz - (dz_corr + tbeta*dr_corr))/ez;
          pulls_dz_all->Fill(pull);
          pulls_dz[icell]->Fill(pull);
        }
        
      }
    }
  }

  RootFile outputFile( outputFilename );
  outputFile.Add( pulls_drphi_all );
  outputFile.Add( pulls_dz_all );
  
  for( auto&& h:pulls_drphi )
  { if( h->GetEntries() ) outputFile.Add( h.release() ); }
    
  for( auto&& h:pulls_dz )
  { if( h->GetEntries() ) outputFile.Add( h.release() ); }

  std::cout << "DistortionsPulls - done" << std::endl;
  
  return outputFilename;
}
