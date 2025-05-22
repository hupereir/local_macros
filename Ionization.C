#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include <TCanvas.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

TH1* nelectrons( TTree* tree, const TString& hname, double gain )
{
      
    // from PHG4MicromegasHitReco
    const double nelectrons_per_gev = 3.73252e+07;
    
    // from Irakli
    const double max_strip_fraction = 0.65;
    const double collection_efficiency = 0.71;
    
    const double electron_charge = 1.6e-19 * 1e15; // fC
    
    // total gain
    const double conversion =  nelectrons_per_gev*gain*max_strip_fraction*collection_efficiency*electron_charge;

    // variable and cut
    const auto var = Form( "DST#EVAL#MicromegasEvaluator_hp::Container._g4hits._eion*%f", conversion );
    TCut cut;
    
    // histogram
    auto h = new TH1F( hname, hname, 100, 0, 500 );
    Utils::TreeToHisto( tree, h->GetName(), var, cut, false );

    auto integrated = Utils::Integrate( h, true, true );
    integrated->GetXaxis()->SetTitle( "q (fc)" );
    integrated->GetYaxis()->SetTitle( "p(q>q_{0})" );
    
    integrated->SetTitle( "" );
    
//     integrated->GetXaxis()->SetTitle( "n_{e-}" );
//     integrated->GetYaxis()->SetTitle( "p(n>n_{e-})" );

    return integrated;
}

TString Ionization(TString tag = TString() )
{

  set_style( false );
  
  if( tag.IsNull() ) tag = "_sHijing_0_12fm";
  const TString inputFile = "DST/CONDOR_Hijing_Micromegas/dst_eval_sHijing_0_12fm/dst_eval*.root";
  const TString pdfFile = Form( "Figures/Ionization%s.pdf", tag.Data() );

  std::cout << "Ionization - inputFile: " << inputFile << std::endl;
  std::cout << "Ionization - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // file manager
  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );

  // primary ionization
  {
    const TString var = "DST#EVAL#MicromegasEvaluator_hp::Container._g4hits._eion*1e6";
    TCut cut;
    
    auto cv = new TCanvas( "cv", "cv", 800, 600 );
    auto h = new TH1F( "h0", "", 100, 0, 15 );
    Utils::TreeToHisto( tree, h->GetName(), var, cut, false );
    
    h->SetTitle("");
    h->GetXaxis()->SetTitle( "E_{ion} (keV)" );
    h->Draw();
    
    gPad->SetLogy( true );
    
    auto mean = h->GetMean();
    auto meanError = h->GetMeanError();
    Draw::VerticalLine( gPad, mean )->Draw();
    
    Draw::PutText( 0.5, 0.75, "MPV = 0.5 keV" );
    Draw::PutText( 0.5, 0.7, Form( "Mean = %.1f keV", mean ) );
    pdfDocument.Add( cv );
    
    auto integrated = Utils::Integrate( h, kTRUE, kTRUE );
    cv =  new TCanvas( "cv", "cv", 800, 600 );
    integrated->GetXaxis()->SetTitle( "E_{ion} (keV)" );
    integrated->GetYaxis()->SetTitle( "p(E>E_{ion})" );
    integrated->SetTitle( "" );
    integrated->Draw( "h" );
    gPad->SetLogy( true );
    Draw::VerticalLine( gPad, mean )->Draw();

    pdfDocument.Add( cv );
  }
  
  // number of electrons per max strip in cluster
  // primary ionization
  {
    auto cv = new TCanvas( "cv", "cv", 800, 600 );

    auto legend = new TLegend( 0.45, 0.8, 0.97, 0.9, "", "NDC" );
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);

    // saturation charge (fC)
    static constexpr double q_sat = 100;

    auto h1 = nelectrons( tree, "h_10000", 10000 );
    const auto val1 = h1->Interpolate( q_sat );
    h1->SetLineColor( 4 );
    h1->Draw( "h" );
    legend->AddEntry( h1, Form( "Gain = 10000, p(q>%.0ffC) = %.2f %%", q_sat, 100.*val1), "L" );
    
    auto h0 = nelectrons( tree, "h_8000", 8000 );
    const auto val0 = h0->Interpolate( q_sat );
    h0->SetLineColor( 2 );
    h0->Draw("same h");
    legend->AddEntry( h0, Form( "Gain = 8000, p(q>%.0ffC) = %.2f %%", q_sat, 100.*val0), "L" );
    
    gPad->SetLogy( true );
    legend->Draw();
    
    // saturation
    Draw::VerticalLine( gPad, q_sat)->Draw();   

    std::cout << "Ionization - sat fraction (G=8000) = " << val0 << std::endl;
    std::cout << "Ionization - sat fraction (G=10000) = " << val1 << std::endl;
    
    pdfDocument.Add( cv );
  }
  
  return pdfFile;
  
}
