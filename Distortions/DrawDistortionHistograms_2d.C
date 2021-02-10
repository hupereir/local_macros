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

static constexpr int isec_rec = 3;
static constexpr double phi_rec = M_PI*isec_rec/6 + M_PI/12;
static constexpr double z_rec = 5;

//_______________________________________________
double Function( double* x, double* par )
{ return par[0] + x[0]*par[1]; }

//_______________________________________________
TString DrawDistortionHistograms_2d()
{

  set_style( false );

  // open reconstructed TFile
  // const TString tag = "_full_realistic_micromegas_mm-empty";
  const TString tag = "_full_realistic_micromegas_truth-empty";
  const auto inputfile = Form( "Rootfiles/Distortions%s.root", tag.Data() );
  std::cout << "DrawDistortionHistograms - inputfile: " << inputfile << std::endl;

  auto f = TFile::Open( inputfile );
  if( !f ) return TString();

  // reconstructed bins
  auto hentries = dynamic_cast<TH3*>(f->Get("hentries_rec"));
  auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hDistortionP_rec"));
  auto hDistortionR_rec = dynamic_cast<TH3*>(f->Get("hDistortionR_rec"));
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

  for( const auto h:{ hentries, hDistortionP_rec, hDistortionR_rec, hDistortionZ_rec } )
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

  // open input map TFile
  const auto mapfile = "distortion_maps/empty.root";
  // const auto mapfile = "distortion_maps/fluct_average-coarse.root";
  // const auto mapfile = "distortion_maps/fluct_average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root";
  // const auto mapfile = "distortion_maps/average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root";
  auto f2 = TFile::Open( mapfile );

  // load histograms
  auto hDRint= dynamic_cast<TH3*>(f2->Get("hIntDistortionR"));
  auto hDPint= dynamic_cast<TH3*>(f2->Get("hIntDistortionP"));
  auto hDZint= dynamic_cast<TH3*>(f2->Get("hIntDistortionZ"));

  // find relevant input bins
  int phibin_in_min = 0;
  int phibin_in_max = 0;
  int zbin_in_min = 0;
  int zbin_in_max = 0;

  for( const auto h: {hDRint, hDPint, hDZint} )
  {
    if( h )
    {
      phibin_in_min = h->GetXaxis()->FindBin( phi_rec_min )+1;
      phibin_in_max = h->GetXaxis()->FindBin( phi_rec_max );

      zbin_in_min = h->GetZaxis()->FindBin( z_rec_min )+1;
      zbin_in_max = h->GetZaxis()->FindBin( z_rec_max );
      break;
    }
  }

  // print bins
  std::cout << "DrawDistortionHistograms - phibin_rec: " << phibin_rec << " phibin_in: (" << phibin_in_min << ", " << phibin_in_max << ")" << std::endl;
  std::cout << "DrawDistortionHistograms - zbin_rec: " << zbin_rec << " zbin_in: (" << zbin_in_min << ", " << zbin_in_max << ")" << std::endl;

  // pdf output
  TString pdfFile( Form( "Figures/DistortionHistograms_2d%s.pdf", tag.Data() ) );
  PdfDocument pdfDocument( pdfFile );

  // load phi histograms
  if( hDistortionP_rec )
  {
    auto cv = new TCanvas( "cv0", "cv0", 900, 900 );
    Draw::DivideCanvas( cv, hDistortionP_rec->GetYaxis()->GetNbins() );

    // project 3D dphi input histogram
    hDPint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
    hDPint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );
    auto projDP = hDPint->Project3D( "y" );
    projDP->Scale( 1./((phibin_in_max-phibin_in_min+1)*(zbin_in_max-zbin_in_min+1)) );

    // project 3D dr histogram
    hDRint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
    hDRint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );
    auto projDR = hDRint->Project3D( "y" );
    projDR->Scale( 1./((phibin_in_max-phibin_in_min+1)*(zbin_in_max-zbin_in_min+1)) );

    // loop over rbins
    for( int ir = 0; ir < hDistortionP_rec->GetYaxis()->GetNbins(); ++ir )
    {

      const auto hname = Form( "residual_2d_drphi_p%i_r%i_z%i", phibin_rec - 1 , ir, zbin_rec - 1 );
      const auto h = dynamic_cast<TH2*>( f->Get( hname ) );
      if( !h ) continue;

      cv->cd( ir+1 );
      h->Draw();

      const auto dphi = hDistortionP_rec->GetBinContent( phibin_rec, ir+1, zbin_rec );
      const auto dr = hDistortionR_rec->GetBinContent( phibin_rec, ir+1, zbin_rec );
      auto fname = Form( "frec_p%i_r%i_z%i", phibin_rec - 1 , ir, zbin_rec - 1 );
      auto f = new TF1( fname, Function, -0.5, 0.5, 2 );
      f->SetParameter(0, dphi );
      f->SetParameter(1, dr );
      f->SetLineColor( 2 );
      f->Draw( "same" );

      // find relevant input r bins
      {
        const auto r_min = hDistortionP_rec->GetYaxis()->GetBinLowEdge(ir+1);
        const auto r_max = hDistortionP_rec->GetYaxis()->GetBinUpEdge(ir+1);
        const auto rbin_min = projDP->GetXaxis()->FindBin( r_min );
        const auto rbin_max = projDP->GetXaxis()->FindBin( r_max );

        float dphi = 0;
        float dr = 0;
        for( int irin = rbin_min; irin <= rbin_max; ++irin )
        {
          dphi += projDP->GetBinContent( irin );
          dr += projDR->GetBinContent( irin );
        }

        dphi /= (rbin_max-rbin_min+1);
        dr /= (rbin_max-rbin_min+1);

        auto fname = Form( "fin_p%i_r%i_z%i", phibin_rec - 1 , ir, zbin_rec - 1 );
        auto f = new TF1( fname, Function, -0.5, 0.5, 2 );
        f->SetParameter(0, dphi );
        f->SetParameter(1, dr );
        f->SetLineColor( 4 );
        f->Draw( "same" );
      }

    }

    pdfDocument.Add(cv);
  }
  
  // get 1D projections with positive and negative alpha separated
  if( true && hDistortionP_rec )
  {
    auto cv = new TCanvas( "cv1", "cv1", 900, 900 );
    Draw::DivideCanvas( cv, hDistortionP_rec->GetYaxis()->GetNbins() );

    // loop over rbins
    for( int ir = 0; ir < hDistortionP_rec->GetYaxis()->GetNbins(); ++ir )
    {

      const auto hname = Form( "residual_2d_drphi_p%i_r%i_z%i", phibin_rec - 1 , ir, zbin_rec - 1 );
      const auto h = dynamic_cast<TH2*>( f->Get( hname ) );
      if( !h ) continue;

      const auto bin = h->GetXaxis()->FindBin( 0. );
      const auto hneg = h->ProjectionY( Form( "%s_neg", hname ), 1, bin-1 );
      const auto hpos = h->ProjectionY( Form( "%s_ps", hname ), bin, h->GetXaxis()->GetNbins() );
      
      cv->cd( ir+1 );
      hneg->SetLineColor(2); 
      // hneg->Scale( 1./hneg->GetEntries() );
      hneg->Draw();
      {
        auto line = Draw::VerticalLine( gPad, hneg->GetMean() );
        line->SetLineColor( 2 );
        line->SetLineWidth(1);
        line->Draw();
      }
      
      hpos->SetLineColor(4);
      // hpos->Scale( 1./hpos->GetEntries() );
      hpos->Draw("same" );
      {
        auto line = Draw::VerticalLine( gPad, hpos->GetMean() );
        line->SetLineColor( 4 );
        line->SetLineWidth(1);
        line->Draw();
      }
    }

    pdfDocument.Add(cv);
  }

  // load z histograms
  if( hDistortionZ_rec )
  {
    auto cv = new TCanvas( "cv2", "cv2", 900, 900 );
    Draw::DivideCanvas( cv, hDistortionZ_rec->GetYaxis()->GetNbins() );

    // project 3D dphi input histogram
    hDZint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
    hDZint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );
    auto projDZ = hDZint->Project3D( "y" );
    projDZ->Scale( 1./((phibin_in_max-phibin_in_min+1)*(zbin_in_max-zbin_in_min+1)) );

    // project 3D dr histogram
    hDRint->GetXaxis()->SetRange( phibin_in_min, phibin_in_max );
    hDRint->GetZaxis()->SetRange( zbin_in_min, zbin_in_max );
    auto projDR = hDRint->Project3D( "y" );
    projDR->Scale( 1./((phibin_in_max-phibin_in_min+1)*(zbin_in_max-zbin_in_min+1)) );

    // loop over rbins
    for( int ir = 0; ir < hDistortionZ_rec->GetYaxis()->GetNbins(); ++ir )
    {

      const auto hname = Form( "residual_2d_dz_p%i_r%i_z%i", phibin_rec - 1 , ir, zbin_rec - 1 );
      const auto h = dynamic_cast<TH2*>( f->Get( hname ) );
      if( !h ) continue;

      cv->cd( ir+1 );
      h->Draw();

      const auto dz = hDistortionZ_rec->GetBinContent( phibin_rec, ir+1, zbin_rec );
      const auto dr = hDistortionR_rec->GetBinContent( phibin_rec, ir+1, zbin_rec );
      auto fname = Form( "frec_dz_p%i_r%i_z%i", phibin_rec - 1 , ir, zbin_rec - 1 );
      auto f = new TF1( fname, Function, -1, 1, 2 );
      f->SetParameter(0, dz );
      f->SetParameter(1, dr );
      f->SetLineColor( 2 );
      f->Draw( "same" );

      // find relevant input r bins
      {
        const auto r_min = hDistortionZ_rec->GetYaxis()->GetBinLowEdge(ir+1);
        const auto r_max = hDistortionZ_rec->GetYaxis()->GetBinUpEdge(ir+1);
        const auto rbin_min = projDZ->GetXaxis()->FindBin( r_min );
        const auto rbin_max = projDZ->GetXaxis()->FindBin( r_max );

        float dz = 0;
        float dr = 0;
        for( int irin = rbin_min; irin <= rbin_max; ++irin )
        {
          dz += projDZ->GetBinContent( irin );
          dr += projDR->GetBinContent( irin );
        }

        dz /= (rbin_max-rbin_min+1);
        dr /= (rbin_max-rbin_min+1);

        auto fname = Form( "fin_dz_p%i_r%i_z%i", phibin_rec - 1 , ir, zbin_rec - 1 );
        auto f = new TF1( fname, Function, -1, 1, 2 );
        f->SetParameter(0, dz );
        f->SetParameter(1, dr );
        f->SetLineColor( 4 );
        f->Draw( "same" );
      }

    }

    pdfDocument.Add(cv);
  }

  return pdfFile;
}
