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

#include "LayerDefines.h"
#include "Fit.C"

//____________________________________________________________________________
namespace
{
  class Result
  {
    public:

    double drphi = 0;
    double dr = 0;
    double drphi_error = 0;
    double dr_error = 0;
    double chisquare = 0;
    double ndf = 0;
  };

  //____________________________________________________________________________
  Result fit_histogram( TH2* h )
  {
    Result out;

    // make a manual fit
    TMatrixD lhs( 2, 2 ); lhs = TMatrixD( TMatrixD::kZero, lhs );
    TMatrixD rhs( 2, 1 ); rhs = TMatrixD( TMatrixD::kZero, rhs );

    const double ybin_width = (h->GetYaxis()->GetXmax() - h->GetYaxis()->GetXmin() )/h->GetNbinsY();

    // loop over histogram bins and fill matrices
    for( Int_t ix = 1; ix <= h->GetNbinsX(); ++ix )
    {
      for( Int_t iy = 1; iy <= h->GetNbinsY(); ++iy )
      {

        const auto x = h->GetXaxis()->GetBinCenter( ix );
        const auto y = h->GetYaxis()->GetBinCenter( iy );
        const auto content = h->GetBinContent( ix, iy );
        if( !content ) continue;

        //uncertainty
        const auto y_error = 1./sqrt(content);
        // const auto y_error = ybin_width/sqrt(content);
        const double erp = square( y_error );

        lhs(0,0) += 1./erp;
        lhs(0,1) += x/erp;
        lhs(1,0) += x/erp;
        lhs(1,1) += square(x)/erp;

        rhs(0,0) += y/erp;
        rhs(1,0) += x*y/erp;

      }
    }

    // invert
    const TMatrixD cov = lhs.InvertFast();
    const TMatrixD result( cov, TMatrixD::kMult, rhs );

    out.drphi = result(0,0);
    out.dr = result(1,0);

    out.drphi_error = std::sqrt( cov(0,0) );
    out.dr_error = std::sqrt( cov(1,1) );

    out.chisquare = 0;
    out.ndf = 0;
    for( Int_t ix = 1; ix <= h->GetNbinsX(); ++ix )
    {
      for( Int_t iy = 1; iy <= h->GetNbinsY(); ++iy )
      {

        const auto x = h->GetXaxis()->GetBinCenter( ix );
        const auto y = h->GetYaxis()->GetBinCenter( iy );
        const auto content = h->GetBinContent( ix, iy );
        if( !content ) continue;

        //uncertainty
        const auto y_error = 1./sqrt(content);
        // const auto y_error = ybin_width/sqrt(content);
        const double erp = square( y_error );

        out.chisquare += square( y - (out.drphi + x*out.dr) )/erp;
        out.ndf += 1;
      }
    }

    return out;
  }

  //____________________________________________________________________________
  Result fit_profile( TProfile* h )
  {
    Result out;

    // make a manual fit
    TMatrixD lhs( 2, 2 ); lhs = TMatrixD( TMatrixD::kZero, lhs );
    TMatrixD rhs( 2, 1 ); rhs = TMatrixD( TMatrixD::kZero, rhs );

    // loop over histogram bins and fill matrices
    for( Int_t ix = 1; ix <= h->GetNbinsX(); ++ix )
    {

      const auto x = h->GetXaxis()->GetBinCenter( ix );
      const auto y = h->GetBinContent( ix );
      const auto y_error = h->GetBinError( ix );

      // skip if error is zero (empty bin)
      if( y_error <= 0 ) continue;

      //uncertainty
      const double erp = square(y_error);

      lhs(0,0) += 1./erp;
      lhs(0,1) += x/erp;
      lhs(1,0) += x/erp;
      lhs(1,1) += square(x)/erp;

      rhs(0,0) += y/erp;
      rhs(1,0) += x*y/erp;

    }

    // invert
    const TMatrixD cov = lhs.InvertFast();
    const TMatrixD result( cov, TMatrixD::kMult, rhs );

    out.drphi = result(0,0);
    out.dr = result(1,0);

    out.drphi_error = std::sqrt( cov(0,0) );
    out.dr_error = std::sqrt( cov(1,1) );

    out.chisquare = 0;
    out.ndf = 0;

    for( Int_t ix = 1; ix <= h->GetNbinsX(); ++ix )
    {

      const auto x = h->GetXaxis()->GetBinCenter( ix );
      const auto y = h->GetBinContent( ix );
      const auto y_error = h->GetBinError( ix );

      // skip if error is zero (empty bin)
      if( y_error <= 0 ) continue;

      //uncertainty
      const double erp = square( y_error);

      out.chisquare += square( y - (out.drphi + x*out.dr) )/erp;
      out.ndf++;
    }

    return out;
  }

}


//____________________________________________________________________________
float delta_phi( float phi )
{
  if( phi >= M_PI ) return phi - 2*M_PI;
  else if( phi < -M_PI ) return phi + 2*M_PI;
  else return phi;
}

//____________________________________________________________________________
double fit_function( double* x, double* par )
{
  const double ta = x[0];
  return par[0]+ta*par[1];
}

//____________________________________________________________________________
TString DeltaRPhi_r_distortions( TString tag = TString() )
{

  set_style( false );
  gStyle->SetPadLeftMargin(0.13);

  constexpr float max_residual = 4;
  constexpr float max_alpha = 0.5;

  // input file
  // if( tag.IsNull() ) tag = "_2k_realistic_full_notpc_noouter" ;
  // if( tag.IsNull() ) tag = "_5k_flat_full_notpc_noouter_notilt" ;
  if( tag.IsNull() ) tag = "_5k_realistic_full_notpc_noouter_notilt_2" ;
  const TString inputFile = Form( "DST/dst_eval%s.root", tag.Data() );

  const TString pdfFile = Form( "Figures/DeltaRPhi_r_distortions%s.pdf", tag.Data() );
  PdfDocument pdfDocument( pdfFile );

  std::cout << "DeltaRPhi2D - inputFile: " << inputFile << std::endl;
  std::cout << "DeltaRPhi2D - pdfFile: " << pdfFile << std::endl;

  constexpr bool fit_profile = false;

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // variable names
  const TString var( "_tracks._clusters._trk_r*delta_phi(_tracks._clusters._trk_phi - _tracks._clusters._phi)" );
  const TString var2d = Form( "%s:-tan( _tracks._clusters._trk_alpha )", var.Data() );
  const TCut momentum_cut;
  const TCut alpha_cut( Form( "fabs( _tracks._clusters._trk_alpha ) < %f", max_alpha ) );
  const TCut residual_cut( Form( "fabs(%s)<%f", var.Data(), max_residual ) );

//   constexpr int nslices = 6;
//   constexpr std::array<int, nslices+1> layers = { 7, 15, 23, 31, 39, 47, 55 };

  constexpr int nslices = 12;
  constexpr std::array<int, nslices+1> layers = { 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55 };

  TCanvas* cv( new TCanvas( "cv", "cv", 800, 800 ) );
  //   Draw::DivideCanvas( cv, nslices, false );
  cv->Divide( 4, 3, 0, 0 );

  // TGraphs
  auto tg = new TGraphErrors();
  auto tg_fit = new TGraphErrors();
  for( int i = 0; i < nslices; ++i )
  {
    const TCut detector_cut( Form( "_tracks._clusters._layer >= %i && _tracks._clusters._layer < %i", layers[i], layers[i+1] ) );

    auto hname = Form( "h%i", i );
    auto h( new TH2F( hname, hname, 100, -max_alpha, max_alpha, 100, -max_residual, max_residual ) );
    Utils::TreeToHisto( tree, hname, var2d, momentum_cut && detector_cut && alpha_cut && residual_cut, false );
    h->SetTitle(TString());
    h->GetXaxis()->SetTitle( "tan(#alpha)" );
    h->GetYaxis()->SetTitle( "r.#Delta#phi_{track-cluster} (cm)" );

    auto pname = Form( "p%i", i );
    auto p( new TProfile( pname, pname, 100, -max_alpha, max_alpha ) );
    Utils::TreeToHisto( tree, pname, var2d, momentum_cut && detector_cut && alpha_cut && residual_cut, false );

    cv->cd( i+1 );
    h->Draw();
    p->SetLineColor(2);
    p->Draw( "same" );

    // calculate mean radius for this slice
    const float r = 0.5*(radius[layers[i]] + radius[layers[i+1]-1]);

    if( false )
    {

      auto fname = Form( "f%i", i );
      auto f = new TF1( fname, fit_function, -max_alpha, max_alpha, 2 );
      f->SetParameter( 0, 0 );

      // define the r transformation (cm)
      static constexpr float deltar_max = 1.5;
      const float dr = -deltar_max*std::cos( 2*M_PI*(r-rmin_tpc)/rlength_tpc );

      // std::cout << "DeltaRPhi_rdistortions - layers: (" << layers[i] << "," << layers[i+1]-1 << ")" << " delta r: " << dr << std::endl;
      f->SetParameter( 1, dr );
      f->SetLineColor(4);
      f->Draw("same");

      tg->SetPoint( i, r, dr );

    } else {

      tg->SetPoint( i, r, 0 );

    }

    // root fit
    auto fname = Form( "f_fit%i", i );
    auto f = new TF1( fname, ::fit_function, -0.5, 0.5, 2 );
    f->SetParameter( 0, 0 );
    f->SetParameter( 1, 0 );

    if( fit_profile ) p->Fit( f, "0QR" );
    else h->Fit( f, "0QR" );
    f->SetLineColor( 2 );
    f->Draw("same");

    const double chisquare_fit = f->GetChisquare();
    const double ndf_fit = f->GetNDF();

    const auto drphi = f->GetParameter(0);
    const auto drphi_error = f->GetParError(0);
    // std::cout << "DeltaRPhi_rdistortions - layers: (" << layers[i] << "," << layers[i+1]-1 << ")" << " delta rphi: " << drphi << " +/- " << drphi_error << std::endl;

    const auto dr = f->GetParameter(1);
    const auto dr_error = f->GetParError(1);
//     std::cout << "DeltaRPhi_rdistortions - layers: (" << layers[i] << "," << layers[i+1]-1 << ")" << " delta r: " << dr << " +/- " << dr_error << std::endl;
//     std::cout << "DeltaRPhi_rdistortions - layers: (" << layers[i] << "," << layers[i+1]-1 << ")" << " chisquare: " << chisquare_fit << " ndf: " << ndf_fit << std::endl;

    tg_fit->SetPoint( i, r, dr );
    tg_fit->SetPointError( i, 0, dr_error );

    // manual histogram fit
    const auto result = fit_profile ? ::fit_profile( p ) : fit_histogram( h );

    std::cout << "DeltaRPhi_rdistortions - layers: (" << layers[i] << "," << layers[i+1]-1 << ")" << " delta r: " << result.dr << " +/- " << result.dr_error << " (calc)" << std::endl;
//     std::cout << "DeltaRPhi_rdistortions - layers: (" << layers[i] << "," << layers[i+1]-1 << ")" << " chisquare: " << result.chisquare << " ndf: " << result.ndf << std::endl;
//     std::cout << "DeltaRPhi_rdistortions - layers: (" << layers[i] << "," << layers[i+1]-1 << ")" << " chisquare ratio: " << chisquare_fit/result.chisquare << std::endl;
//     std::cout << std::endl;

  }

  pdfDocument.Add( cv );

  // draw tgraphs
  cv = new TCanvas( "cv1", "cv1", 800, 800 );
  tg->SetMarkerStyle( 20 );
  tg->SetMarkerColor( 4 );
  tg->SetLineColor( 4 );
  tg->GetXaxis()->SetTitle( "r (cm)" );
  tg->GetYaxis()->SetTitle( "#Deltar (cm)" );
  tg->SetMinimum( -2.5 );
  tg->SetMaximum( 2.5 );
  tg->Draw("AP");

  tg_fit->SetMarkerStyle( 20 );
  tg_fit->SetMarkerColor( 2 );
  tg_fit->SetLineColor( 2 );
  tg_fit->Draw("P");
  cv->Update();

  auto legend = new TLegend( 0.16, 0.82, 0.60, 0.93, "", "NDC" );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->Draw();
  legend->AddEntry( tg, "input distortion", "AP" );
  legend->AddEntry( tg_fit, "reconstructed distortion", "AP" );

  pdfDocument.Add( cv );

  return pdfFile;

}
