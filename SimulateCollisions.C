#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Stream.h>
#include <RootUtil/Utils.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h> 

R__LOAD_LIBRARY(libgslcblas.so)
R__LOAD_LIBRARY(libgsl.so)

void SimulateCollisions()
{
  
  set_style( false );
  
  PdfDocument pdfDocument( "Figures/timestamps.pdf" );
  
  // random generator
  auto rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set( rng, 1 );

  // time interval (s) between two bunch crossing
  /* value copied from generators/phhepmc/Fun4AllHepMCPileupInputManager.cc */
  static constexpr double deltat_crossing = 106e-9; 
  std::cout << "SimulateCollisions - deltat_crossing: " << deltat_crossing << std::endl;
  
  // triggered rate (Hz)
  static constexpr double trigger_rate = 1e5;
  
  // mean number of collision per crossing
  static constexpr double mu = trigger_rate*deltat_crossing;
  
  // print configuration
  std::cout << "SimulateCollisions - deltat_crossing: " << deltat_crossing << std::endl;
  std::cout << "SimulateCollisions - trigger_rate: " << trigger_rate << std::endl;
  std::cout << "SimulateCollisions - mu: " << mu << std::endl;
  
  // number of requested triggers (should correspond approximately to 1/10 second of data taking)
  static constexpr uint ntrigtot = 1e4;

  // put interval between collisions into a histogram
  auto h_deltat = new TH1F( "deltat", "", 10./(trigger_rate*deltat_crossing), 0, 10./trigger_rate );
  auto h_prob = new TH1F( "prob", "", 5, 0, 5 );
  
  // trigger time version collision number
  using trigger_time_t = std::tuple<int64_t,double>;
  std::array<trigger_time_t,ntrigtot> trigger_times;
    
  // running collision time
  int64_t bunchcrossing = 0;
  double time = 0;

  // keep track of the last collision time
  double previous_trigger_time = 0;
  
  // generate triggers
  for( int itrig = 0; itrig < ntrigtot; )
  {
    
    ++ bunchcrossing;
    time += deltat_crossing;
    
    auto ntrig = gsl_ran_poisson( rng, mu );
    h_prob->Fill( ntrig );
    
    for( uint i = 0; i < ntrig; ++i )
    { 
      std::cout << "SimulateCollisions - trigger number: " << itrig << " bunch crossing: " << bunchcrossing << " time: " << time << std::endl;
      trigger_times[itrig++] = {bunchcrossing, time};      
      h_deltat->Fill( time - previous_trigger_time );
      previous_trigger_time = time;
    }        
  }
  
  // printout last trigger time
  std::cout << "SimulateCollisions - last trigger time: " << time << std::endl;
  
  {
    // draw deltat histogram
    auto cv = new TCanvas( "cv", "cv", 800, 800 );
    cv->SetRightMargin( .12 );
    h_deltat->GetXaxis()->SetTitle( "#Delta t (s)" );
    h_deltat->GetXaxis()->SetMaxDigits( 1 );
    h_deltat->Draw();
    gPad->SetLogy( true );
    Draw::PutText( 0.4, 0.8, Form( "Mean = %.1g #pm %.1g", h_deltat->GetMean(), h_deltat->GetMeanError() ) );
    
    // mean interval (should be 1/trigger_rate)
    pdfDocument.Add( cv );
    std::cout << "SimulateCollisions - mean interval between collisions: " << h_deltat->GetMean() << std::endl;
  }
  
  {
    // draw probability histogram
    auto cv = new TCanvas( "cv1", "cv1", 800, 800 );
    // cv->SetRightMargin( .12 );
    h_prob->GetXaxis()->SetTitle( "trigger/bunch cross" );
    h_prob->Draw();
    gPad->SetLogy( true );    

    Draw::PutText( 0.4, 0.8, Form( "Mean = %.3g #pm %.3g", h_prob->GetMean(), h_prob->GetMeanError() ) );
    
    // fraction of bunch crossing with more than 1 collision
    pdfDocument.Add( cv );
    std::cout << "SimulateCollisions - in-bunch pileup: " << 1. - (h_prob->GetBinContent(1) + h_prob->GetBinContent(2))/h_prob->GetEntries() << std::endl;
  }
      
  gsl_rng_free(rng);
  
  // write all timestamps to a file
  std::ofstream out( "timestamps.txt" );
  out << "// bunchcrossin id; time (ns)" << std::endl;
  out << "// assuming 106ns between bunches" << std::endl;
  for( const auto& [bunchcrossing,time]:trigger_times ) { out << bunchcrossing << " " << time*1e9 << std::endl; }
  out.close();

}
