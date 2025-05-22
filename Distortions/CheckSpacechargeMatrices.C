#include <RootUtil/Draw.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>
#include <RootUtil/FileManager.h>
#include <tpccalib/TpcSpaceChargeMatrixContainerv1.h>
#include <frog/FROG.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libtpccalib.so)

namespace
{
  // phi range
  static constexpr float m_phimin = 0;
  static constexpr float m_phimax = 2.*M_PI;

  // TODO: could try to get the r and z range from TPC geometry
  // r range
  static constexpr float m_rmin = 20;
  static constexpr float m_rmax = 78;

  // z range
  static constexpr float m_zmin = -105.5;
  static constexpr float m_zmax = 105.5;
  
  
  static constexpr double z_rec = 5;

  double get_sector_phi( int isec ) { return isec*M_PI/6 + M_PI/12; } 
}

//_______________________________________________
TString CheckSpacechargeMatrices()
{

  const std::string inputfile = "tpcspacechargematrices.list";

  // create container
  
  auto container = new TpcSpaceChargeMatrixContainerv1;
  
  FROG frog;
  std::ifstream in( inputfile.c_str() );
  const int max_files = 1000;
  int ifile = 0;
  while( !(in.rdstate()&std::ios::failbit ) )
  {
    if( max_files > 0 && ++ifile >= max_files ) break;
    
    std::string line;
    std::getline( in, line );
    if( line.empty() ) continue;
    if( line.substr(0,2) == "//" )
    {
      std::cout << "CheckSpacechargeMatrices - skipping " << line << std::endl;
      continue;
    }

    const auto long_filename = frog.location( line );
    std::cout << "CheckSpacechargeMatrices - processing: " << long_filename << std::endl;
    
    // save everything to root file
    std::unique_ptr<TFile> inputfile( TFile::Open( long_filename ) );
    std::unique_ptr<TpcSpaceChargeMatrixContainer> source( dynamic_cast<TpcSpaceChargeMatrixContainer*>( inputfile->Get( "TpcSpaceChargeMatrixContainer" ) ) );
    container->add( *source );
  }
  
  // pdf output
  const TString pdfFile = "Figures/TpcSpaceChargeMatrices.pdf";
  const TString outputfile = "Rootfiles/TpcSpaceChargeMatrices.root";
  PdfDocument pdfDocument( pdfFile );

  /
  
  container->identify();
  
  // get grid dimensions from matrix container
  int phibins = 0;
  int rbins = 0;
  int zbins = 0;
  container->get_grid_dimensions( phibins, rbins, zbins );
  
  // fill entries
  auto hentries( new TH3F( "hentries_rec", "hentries_rec", phibins, m_phimin, m_phimax, rbins, m_rmin, m_rmax, zbins, m_zmin, m_zmax ) );

  // set axis labels
  for( const auto& h:{ hentries } )
  {
    h->GetXaxis()->SetTitle( "#phi (rad)" );
    h->GetYaxis()->SetTitle( "r (cm)" );
    h->GetZaxis()->SetTitle( "z (cm)" );
  }
  
  // loop over bins
  for( int iphi = 0; iphi < phibins; ++iphi )
    for( int ir = 0; ir < rbins; ++ir )
    for( int iz = 0; iz < zbins; ++iz )
  {
    // get cell index
    const auto icell = container->get_cell_index( iphi, ir, iz );
    const auto cell_entries = container->get_entries(icell);
    hentries->SetBinContent( iphi+1, ir+1, iz+1, cell_entries );
  }
  
  {
    // entries in r, phi plane
    TCanvas* cv( new TCanvas( "cv1", "cv1", 800, 800 ) );
    cv->SetRightMargin( 0.24 );
    
    auto zbin_rec = hentries->GetZaxis()->FindBin( z_rec );
    hentries->GetZaxis()->SetRange( zbin_rec, zbin_rec ); // z axis
    auto proj = hentries->Project3D( "yx" );
    proj->SetTitle("");
    proj->GetZaxis()->SetTitle( "entries" );
    proj->GetZaxis()->SetTitleOffset( 2.1 );
    proj->Draw( "colz" );
    
    // draw sector boundaries
    for( int i = 0; i < 12; ++i )
    {
      const auto phi( get_sector_phi(i) + M_PI/12 );
      Draw::VerticalLine( cv, phi )->Draw();
    }
    
    {
      const int isec = 3;
      const auto phi_rec = get_sector_phi(isec);
      auto line = Draw::VerticalLine( cv, phi_rec );
      line->SetLineColor( 2 );
      line->Draw();
    }
    
    pdfDocument.Add(cv);
    
  }
  
  // also agreggate files to single output
  
  // save matrix container in output file
  if( m_matrix_container )
  {
    std::unique_ptr<TFile> outputfile( TFile::Open( outputfile, "RECREATE" ) );
    outputfile->cd();
    m_matrix_container->Write( "TpcSpaceChargeMatrixContainer" );
  }
  
  return pdfFile; 
}
