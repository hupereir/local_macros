#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/RootFile.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>

const std::string status = "Preliminary";

TGraphErrors* get_ref_efficiency( double eff_min = 0, double eff_max = 1)
{

  // reference efficiency from Bob's lab
  const double hv[] = { 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440 };
  const double eff[] = { 0.09, 0.178, 0.294, .45, 0.603, 0.77, 0.864, 0.94, 0.97, 0.98, 0.99 };
  const unsigned int npoints = 11;
  
  auto tg = new TGraphErrors();
  tg->SetTitle("");
  tg->GetXaxis()->SetTitle( "HV_{resist}" );
  tg->GetYaxis()->SetTitle( "eff" );
  tg->SetMarkerStyle( 20 );
  tg->SetMarkerColor( 2 );

  for( int i = 0; i < npoints; ++i )
  {
    double eff_corrected = eff_min*(1.-eff[i]) + eff_max*eff[i];
    tg->SetPoint(i, hv[i], eff_corrected );
  }
  
  return tg;
}

void DrawRawDataEfficiency()
{
  #if true
  std::vector<TString> inputfiles = {
//     "Rootfiles/RawDataEfficiency_phi-D400_4sigma.root",
//     "Rootfiles/RawDataEfficiency_phi-D400_4sigma_charge_cut.root",
//     "Rootfiles/RawDataEfficiency_phi-D400_4sigma_size_cut.root",
//     "Rootfiles/RawDataEfficiency_phi-D400_4sigma_charge_size_cut.root",
//     "Rootfiles/RawDataEfficiency_phi-D400_4.5sigma.root",
//     "Rootfiles/RawDataEfficiency_phi-D400_4.5sigma_charge_cut.root",
//     "Rootfiles/RawDataEfficiency_phi-D400_4.5sigma_size_cut.root",
//     "Rootfiles/RawDataEfficiency_phi-D400_4.5sigma_charge_size_cut.root",
    "Rootfiles/RawDataEfficiency_phi-D400_nocut.root",
    "Rootfiles/RawDataEfficiency_phi-D400_charge_cut.root",
    "Rootfiles/RawDataEfficiency_phi-D400_size_cut.root",
    "Rootfiles/RawDataEfficiency_phi-D400_charge_size_cut.root"
  };
  
  const std::vector<TString> detectors = 
  {
    "SCIP",
    "NCOP",
    "SEIP",
    "NEIP", 
    "SWIP",
    "NWIP"
  };
  
  const bool is_phi = true;
  #else
  std::vector<TString> inputfiles = {
//     "Rootfiles/RawDataEfficiency_z-D400_4sigma.root",
//     "Rootfiles/RawDataEfficiency_z-D400_4sigma_charge_cut.root",
//     "Rootfiles/RawDataEfficiency_z-D400_4sigma_size_cut.root",
//     "Rootfiles/RawDataEfficiency_z-D400_4sigma_charge_size_cut.root",
//     "Rootfiles/RawDataEfficiency_z-D400_4.5sigma.root",
//     "Rootfiles/RawDataEfficiency_z-D400_4.5sigma_charge_cut.root",
//     "Rootfiles/RawDataEfficiency_z-D400_4.5sigma_size_cut.root",
//     "Rootfiles/RawDataEfficiency_z-D400_4.5sigma_charge_size_cut.root",
    "Rootfiles/RawDataEfficiency_z-D400_nocut.root",
    "Rootfiles/RawDataEfficiency_z-D400_charge_cut.root",
    "Rootfiles/RawDataEfficiency_z-D400_size_cut.root",
    "Rootfiles/RawDataEfficiency_z-D400_charge_size_cut.root"
  };
  const TString detector = "SCIZ"; 

  const std::vector<TString> detectors = 
  {
    "SCIZ",
    "NCOZ",
    "SEIZ",
    "NEIZ", 
    "SWIZ",
    "NWIZ"
  };
  const bool is_phi = false;
  #endif
  
  for( const auto& detector:detectors )
  {
    
    std::vector<double> hv(12, 0);
    std::vector<double> value_mean( 12, 0 );
    std::vector<double> error_mean( 12, 0 );
    std::vector<double> value_min( 12, 0 );
    std::vector<double> value_max( 12, 0 );
    
    // create plot
    TCanvas* cv0 = new TCanvas( "cv0", "cv", 980, 900 );
    {
      auto h = new TH1I( "h", "", 100, 310, 440 );
      h->SetMinimum(0);
      h->SetMaximum(1);
      h->GetXaxis()->SetTitle( "resist HV (V)" );
      h->GetYaxis()->SetTitle( "efficiency" );
      h->Draw();
      gPad->SetTopMargin( 0.07 );
      gPad->SetLeftMargin( 0.14);
      gPad->Update();
    }
    
    // fill values from files
    bool first = true;
    int color = 1;
    for( const auto& inputfile:inputfiles )
    {
      std::cout << "DrawRawDataEfficiency - file: " << inputfile << std::endl;
      auto tfile = std::make_unique<TFile>( inputfile, "READ" );
      auto tg = static_cast<TGraphErrors*>( tfile->Get( Form( "tge_%s", detector.Data() ) ) );
      
      if( first ) color = tg->GetMarkerColor();
      
      for( int i =0; i < tg->GetN(); ++i )
      {
        const double x = tg->GetPointX(i);
        const double y = tg->GetPointY(i);
        const double ey = tg->GetErrorY(i);
        std::cout << "x: " << x << " y: " << y << " ey: " << ey << std::endl;
        
        // save/check x axis
        if( hv[i] == 0 ) hv[i] = x;
        else assert( hv[i]==x );
        
        // update min/max
        if( first || y < value_min[i] ) value_min[i] = y;
        if( first || y > value_max[i] ) value_max[i] = y;
        
        // update mean and error
        value_mean[i] += y;
        error_mean[i] += ey;
      }
      
      first = false;

      tg->SetMarkerSize(2);
      tg->Draw("P");
      
    }
    
    for( auto&& value:value_mean ) { value /= inputfiles.size(); }
    for( auto&& error:error_mean ) { error /= inputfiles.size(); }

    cv0->SaveAs( Form( "Figures/RawDataEfficiency_%s_raw.pdf", detector.Data() ) );
    cv0->SaveAs( Form( "Figures/RawDataEfficiency_%s_raw.png", detector.Data() ) );
   
    // create plot
    TCanvas* cv = new TCanvas( "cv", "cv", 980, 900 );
    
    auto h = new TH1I( "h", "", 100, 310, 440 );
    h->SetMinimum(0);
    h->SetMaximum(1);
    h->GetXaxis()->SetTitle( "resist HV (V)" );
    h->GetYaxis()->SetTitle( "efficiency" );
    h->Draw();
    gPad->SetTopMargin( 0.07 );
    gPad->SetLeftMargin( 0.14);
    
    gPad->Update();
    
    auto tge = new TGraphErrors();
    tge->SetMarkerStyle(20);
    tge->SetMarkerColor(color);
    tge->SetLineColor( color );
    tge->SetMarkerSize(2);    
    double eff_min = 1.;
    double eff_max = 0;
    for( size_t i = 0; i < hv.size(); ++i )
    { 
      eff_min = std::min( eff_min, value_mean[i] );
      eff_max = std::max( eff_max, value_mean[i] );

      tge->SetPoint( i,  hv[i], value_mean[i] );
      tge->SetPointError( i, 0, error_mean[i] );
    }
    
    tge->Draw("P");
    
    auto tg_ref = get_ref_efficiency(eff_min, eff_max);
    tg_ref->Draw("P");
    
    
    for( size_t i = 0; i < hv.size(); ++i )
    {
      Draw::SetBoxFillStyle(3002);
      Draw::DrawBox( hv[i], value_mean[i], 2, (value_max[i]-value_min[i] )/std::sqrt(12), color );
    }
    
    {
      auto text = new TLatex;
      text->SetNDC( true );
      text->SetTextSize( 0.045 );
      text->DrawLatex( 0.77, 0.95, "#it{08/30/2023}" );
    }
    
    {
      auto text = new TPaveText(0.17,0.73,0.53,0.88, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText( Form( "#it{#bf{sPHENIX}} %s", status.c_str() ));
      text->AddText("Au+Au #sqrt{s_{NN}} = 200 GeV");
      text->AddText(Form( "TPOT %s (%s strips)", detector.Data(), is_phi ? "#Phi":"#it{z}" ) );
      text->Draw();
    }
    
    cv->SaveAs( Form( "Figures/RawDataEfficiency_%s.pdf", detector.Data() ) );
    cv->SaveAs( Form( "Figures/RawDataEfficiency_%s.png", detector.Data() ) );
  }
}
