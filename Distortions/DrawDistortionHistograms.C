#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

namespace
{
  int GetRangeBins( TAxis* axis ) 
  { return axis->GetLast() - axis->GetFirst() + 1; }
 
  constexpr double get_sector_phi( int isec ) 
  { return isec*M_PI/6 + M_PI/12; }
 
}

// specify bins for which residuals are retreived
static constexpr double phi_rec =  get_sector_phi(3);
static constexpr double z_rec = 5.;

//_______________________________________________
TString DrawDistortionHistograms()
{

  set_style( false );

  // open TFile
  const TString tag = "_full_realistic_micromegas_all-coarse";
  const auto inputfile = Form( "Rootfiles/Distortions%s.root", tag.Data() );

  std::cout << "DrawDistortionHistograms - inputfile: " << inputfile << std::endl;
  auto f = TFile::Open( inputfile );
  if( !f ) return TString();

  // configuration
  const bool use_range = false;
  
  // load histograms
  auto hentries = dynamic_cast<TH3*>(f->Get("hentries_rec"));
  auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hDistortionP_rec"));
  auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hDistortionZ_rec"));

  // find relevant reconstructed bins
  std::cout << "DrawDistortionHistograms - phi_rec: " << phi_rec << std::endl;
  std::cout << "DrawDistortionHistograms - z_rec: " << z_rec << std::endl;

  int phibin_rec = 0;
  double phi_rec_min = phi_rec;
  double phi_rec_max = phi_rec;

  int zbin_rec = 0;
  double z_rec_min = z_rec;
  double z_rec_max = z_rec;

  for( const auto h:{ hentries, hDistortionP_rec, hDistortionZ_rec } )
  {
    if( h )
    {
      phibin_rec = h->GetXaxis()->FindBin( phi_rec );
      phi_rec_min = h->GetXaxis()->GetBinLowEdge( phibin_rec );
      phi_rec_max = h->GetXaxis()->GetBinUpEdge( phibin_rec );

      zbin_rec = h->GetZaxis()->FindBin( z_rec );
      z_rec_min = h->GetZaxis()->GetBinLowEdge( zbin_rec );
      z_rec_max = h->GetZaxis()->GetBinUpEdge( zbin_rec );
      break;
    }
  }

  // const auto mapfile = "distortion_maps/empty.root";
  const auto mapfile = "distortion_maps/fluct_average-coarse.root";
  auto f2 = TFile::Open( mapfile );
  if( !f2 ) return TString();

  // load histograms
  auto hDPint= dynamic_cast<TH3*>(f2->Get("hIntDistortionP"));
  auto hDZint= dynamic_cast<TH3*>(f2->Get("hIntDistortionZ"));

  // find relevant input bins
  int phibin_in = 0;
  int phibin_in_min = 0;
  int phibin_in_max = 0;
  
  int zbin_in = 0;
  int zbin_in_min = 0;
  int zbin_in_max = 0;

  for( const auto h: {hDPint, hDZint} )
  {
    if( h )
    {
      phibin_in = h->GetXaxis()->FindBin( phi_rec );
      phibin_in_min = h->GetXaxis()->FindBin( phi_rec_min )+1;
      phibin_in_max = h->GetXaxis()->FindBin( phi_rec_max );

      zbin_in = h->GetXaxis()->FindBin( z_rec );
      zbin_in_min = h->GetZaxis()->FindBin( z_rec_min )+1;
      zbin_in_max = h->GetZaxis()->FindBin( z_rec_max );
      break;
    }
  }

  // print bins
  std::cout << "DrawDistortionHistograms - phibin_rec: " << phibin_rec << " phibin_in: (" << phibin_in_min << ", " << phibin_in_max << ")" << std::endl;
  std::cout << "DrawDistortionHistograms - zbin_rec: " << zbin_rec << " zbin_in: (" << zbin_in_min << ", " << zbin_in_max << ")" << std::endl;

  // pdf output
  TString pdfFile( Form( "Figures/DistortionHistograms%s.pdf", tag.Data() ) );
  PdfDocument pdfDocument( pdfFile );

  using histogram_struct_t = std::tuple<TH3*, TH3*, TString>;
  std::array<histogram_struct_t, 2> histogram_structs = {{
    { hDistortionP_rec, hDPint, "residual_drphi_p%i_r%i_z%i" },
    { hDistortionZ_rec, hDZint, "residual_dz_p%i_r%i_z%i" }
  }};

  int icv = 0;
  for( const auto& [hdistortion_rec, hint, hname]:histogram_structs )
  {  
    if( !hdistortion_rec ) continue;
    if( !hint ) continue;
    
    // create canvas
    const auto cvname = Form( "cv_%i", ++icv );
    auto cv = new TCanvas( cvname, cvname, 900, 900 );
    
    // divide
    const auto nrbins = hdistortion_rec->GetYaxis()->GetNbins();
    Draw::DivideCanvas( cv, nrbins );

    // project 3D input histogram
    if( use_range )
    {
      hint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
      hint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );
    } else {
      hint->GetXaxis()->SetRange( phibin_in, phibin_in );
      hint->GetZaxis()->SetRange( zbin_in, zbin_in );
    }
    
    auto proj = hint->Project3D( "y" );
    proj->Scale( 1./(GetRangeBins(hint->GetXaxis())*GetRangeBins(hint->GetZaxis())) );

    // loop over rbins
    for( int ir = 0; ir < nrbins; ++ir )
    {
      
      // draw residual distribution
      const auto fullname = Form( hname.Data(), phibin_rec - 1 , ir, zbin_rec - 1 );
      const auto h = dynamic_cast<TH1*>( f->Get( fullname ) );
      if( !h ) continue;

      cv->cd( ir+1 );
      h->Draw();

      // draw reconstructed distortion
      const auto value = hdistortion_rec->GetBinContent( phibin_rec, ir+1, zbin_rec );
      auto line = Draw::VerticalLine( gPad, value );
      line->SetLineColor( 2 );
      line->Draw();

      // draw input distortion
      {
        float value_in = 0;
        if( use_range )
        {
          
          const auto r_min = hdistortion_rec->GetYaxis()->GetBinLowEdge(ir+1);
          const auto r_max = hdistortion_rec->GetYaxis()->GetBinUpEdge(ir+1);
          const auto rbin_min = proj->GetXaxis()->FindBin( r_min );
          const auto rbin_max = proj->GetXaxis()->FindBin( r_max );
          
          for( int irin = rbin_min; irin <= rbin_max; ++irin )
          { value_in += proj->GetBinContent( irin ); }
        
          value_in /= (rbin_max-rbin_min+1);
          
        } else {
          
          const auto r = hdistortion_rec->GetYaxis()->GetBinCenter(ir+1);
          const auto irin = proj->GetXaxis()->FindBin( r );
          value_in = proj->GetBinContent( irin );

        }

        auto line = Draw::VerticalLine( gPad, value_in );
        line->SetLineColor( 4 );
        line->Draw();
      }
    }

    pdfDocument.Add(cv);

  }

  return pdfFile;
}
