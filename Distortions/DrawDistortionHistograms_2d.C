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
  { return isec*M_PI/6; }
 
}

// specify bins for which residuals are retreived
static constexpr int isec = 9;
static constexpr double phi_rec =  get_sector_phi(isec);
static constexpr double z_rec = 5.;

//_______________________________________________
double Function( double* x, double* par )
{ return par[0] + x[0]*par[1]; }

//_______________________________________________
TString DrawDistortionHistograms_2d()
{

  set_style( false );

  // open TFile 
//   const TString tag = "_full_realistic_micromegas_mm_genfit_nodistortion";
//   const auto inputmap = Form( "Rootfiles/Distortions%s.root", tag.Data() );
//   const auto inputhistogram = "DST/TpcResiduals_genfit_truth_no_distortion.root";

//   // open TFile 
//   const TString tag = "_full_realistic_micromegas_mm_genfit_nodistortion";
//   const auto inputmap = Form( "Rootfiles/Distortions%s.root", tag.Data() );
//   const auto inputhistogram = "DST/CONDOR_realistic_micromegas/dst_reco_genfit_truth_notpc_nodistortion/TpcResiduals_*.root";

  const TString tag = "_full_realistic_micromegas_mm_acts_nodistortion-test3";
  // const TString tag = "_full_realistic_micromegas_mm_acts_nodistortion";
  const auto inputmap = Form( "Rootfiles/Distortions%s.root", tag.Data() );
  const auto inputhistogram = "DST/CONDOR_realistic_micromegas/dst_reco_acts_truth_notpc_nodistortion-test3/TpcResiduals_*.root";

  std::cout << "DrawDistortionHistograms - inputmap: " << inputmap << std::endl;
  std::cout << "DrawDistortionHistograms - inputhistogram: " << inputhistogram << std::endl;
  FileManager filemanager( inputmap );
  FileManager histogrammanager( inputhistogram );

  // configuration
  const bool use_range = false;
  const bool do_fit = true;
  
  // reconstructed bins
  auto hentries = dynamic_cast<TH3*>(filemanager.GetHistogram("hentries_rec"));
  auto hDistortionP_rec = dynamic_cast<TH3*>(filemanager.GetHistogram("hDistortionP_rec"));
  auto hDistortionR_rec = dynamic_cast<TH3*>(filemanager.GetHistogram("hDistortionR_rec"));
  auto hDistortionZ_rec = dynamic_cast<TH3*>(filemanager.GetHistogram("hDistortionZ_rec"));
  
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
  // const auto mapfile = "distortion_maps-new/average_minus_static_distortion_coarse.root";
  const auto mapfile = "distortion_maps-new/empty.root";
  auto f2 = TFile::Open( mapfile );
  if( !f2 ) return TString();

  // load histograms
  auto hDRint= dynamic_cast<TH3*>(f2->Get("hIntDistortionR_posz"));
  auto hDPint= dynamic_cast<TH3*>(f2->Get("hIntDistortionP_posz"));
  auto hDZint= dynamic_cast<TH3*>(f2->Get("hIntDistortionZ_posz"));

  // find relevant input bins
  int phibin_in = 0;
  int phibin_in_min = 0;
  int phibin_in_max = 0;
  
  int zbin_in = 0;
  int zbin_in_min = 0;
  int zbin_in_max = 0;

  for( const auto h: {hDPint, hDRint, hDZint} )
  {
    if( h )
    {
      phibin_in = h->GetXaxis()->FindBin( phi_rec );
      phibin_in_min = h->GetXaxis()->FindBin( phi_rec_min )+1;
      phibin_in_max = h->GetXaxis()->FindBin( phi_rec_max );

      zbin_in = h->GetZaxis()->FindBin( z_rec );
      zbin_in_min = h->GetZaxis()->FindBin( z_rec_min )+1;
      zbin_in_max = h->GetZaxis()->FindBin( z_rec_max );
      break;
    }
  }
  
  // print bins
  std::cout << "DrawDistortionHistograms - phibin_rec: " << phibin_rec << " phibin_in: (" << phibin_in << "," << phibin_in_min << ", " << phibin_in_max << ")" << std::endl;
  std::cout << "DrawDistortionHistograms - zbin_rec: " << zbin_rec << " zbin_in: (" << zbin_in << "," << zbin_in_min << ", " << zbin_in_max << ")" << std::endl;

  // pdf output
  TString pdfFile( Form( "Figures/DistortionHistograms_2d%s_%i.pdf", tag.Data(), isec ) );
  PdfDocument pdfDocument( pdfFile );
  
  using histogram_struct_t = std::tuple<TH3*, TH3*, TString>;
  std::array<histogram_struct_t, 2> histogram_structs = {{
    { hDistortionP_rec, hDPint, "residual_2d_drphi_p%i_r%i_z%i" },
    { hDistortionZ_rec, hDZint, "residual_2d_dz_p%i_r%i_z%i" }
  }};


  int icv = 0;
  for( const auto& [hdistortion_rec, hint, hname]:histogram_structs )
  {
    
    std::cout << "DrawDistortionHistograms_2d - hname: " << hname << std::endl;
    
    if( !hdistortion_rec ) continue;
    if( !hint ) continue;
    
    // create canvas
    const auto cvname = Form( "cv_%i", ++icv );
    auto cv = new TCanvas( cvname, cvname, 900, 900 );
    
    // divide
    const auto nrbins = hdistortion_rec->GetYaxis()->GetNbins();
    Draw::DivideCanvas( cv, nrbins );

    bool first = true;
    
    // project 3D histograms into relevant 1D histogram vs r
    auto get_projection = [&]( TH3* hint ) 
    {
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
      return proj;
    };
    
    auto proj = get_projection( hint );
    auto projDR = get_projection( hDRint );

    // loop over rbins
    TLegend* legend = nullptr;
    for( int ir = 0; ir < nrbins; ++ir )
    {
      std::cout << "DrawDistortionHistograms_2d - ir: " << ir << std::endl;
     
      // draw residual distribution
      const auto fullname = Form( hname.Data(), phibin_rec - 1 , ir, zbin_rec - 1 );
      const auto h = dynamic_cast<TH2*>( histogrammanager.GetHistogram( fullname ) );
      if( !h ) continue;

      cv->cd( ir+1 );
      h->Draw( "col" );

      // draw legend
      if( first )
      {
        legend = new TLegend( 0.16, 0.15, 0.60, 0.3, "", "NDC" );
        legend->SetFillColor(0);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->Draw();
      }
      
      // draw input distortion
      {
        float value_in = 0;
        float dr_in = 0;
        if( use_range )
        {
          
          const auto r_min = hdistortion_rec->GetYaxis()->GetBinLowEdge(ir+1);
          const auto r_max = hdistortion_rec->GetYaxis()->GetBinUpEdge(ir+1);
          const auto rbin_min = proj->GetXaxis()->FindBin( r_min );
          const auto rbin_max = proj->GetXaxis()->FindBin( r_max );
          
          for( int irin = rbin_min; irin <= rbin_max; ++irin )
          {
            value_in += proj->GetBinContent( irin );
            dr_in += projDR->GetBinContent( irin );
          }
        
          value_in /= (rbin_max-rbin_min+1);
          dr_in /= (rbin_max-rbin_min+1);
          
        } else {
          
          const auto r = hdistortion_rec->GetYaxis()->GetBinCenter(ir+1);
          const auto irin = proj->GetXaxis()->FindBin( r );
          value_in = proj->GetBinContent( irin );
          dr_in = projDR->GetBinContent( irin );

        }

        auto fname = Form( "fin_p%i_r%i_z%i", phibin_rec - 1 , ir, zbin_rec - 1 );
        auto f = new TF1( fname, Function, -0.5, 0.5, 2 );
        f->SetParameter(0, value_in );
        f->SetParameter(1, dr_in );
        f->SetLineColor( 4 );
        f->Draw( "same" );
        
        if( first ) legend->AddEntry( f, "input", "L" );
      
      }
      
      {
        // draw reconstructed distortion
        const auto value = hdistortion_rec->GetBinContent( phibin_rec, ir+1, zbin_rec );
        const auto dr = hDistortionR_rec->GetBinContent( phibin_rec, ir+1, zbin_rec );
        const auto fname = Form( "frec_p%i_r%i_z%i", phibin_rec - 1 , ir, zbin_rec - 1 );
        auto f = new TF1( fname, Function, -0.5, 0.5, 2 );
        f->SetParameter(0, value );
        f->SetParameter(1, dr );
        f->SetLineColor( 2 );
        f->Draw( "same" );

        if( first ) legend->AddEntry( f, "reconstructed", "L" );
      }
      
      if( do_fit )
      {
        // fit
        const auto fname = Form( "ffit_p%i_r%i_z%i", phibin_rec - 1 , ir, zbin_rec - 1 );
        auto f = new TF1( fname, Function, -0.5, 0.5, 2 );
        f->SetParameter(0,0);
        f->SetParameter(1,0);
        h->Fit( f, "0QR" );
        f->SetLineColor( 6 );
        f->Draw( "same" );

        if( first ) legend->AddEntry( f, "fit", "L" );
      }
        
      first = false;
    }

    pdfDocument.Add(cv);
  }

  return pdfFile;
}
