#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>

#include <TPRegexp.h>
#include <TString.h>

#include <fstream>
#include <string>
#include <map>
#include <vector>

R__LOAD_LIBRARY(libRootUtilBase.so)

namespace
{
  using TimeArray = std::vector<double>;
  using ModuleMap = std::map<std::string, TimeArray>;
}

//________________________________________________________________________________
void ParseFile( const std::string& filename, ModuleMap& module_times )
{
  std::cout << "ParseFile - filename: " << filename << std::endl;
  std::ifstream in( filename.c_str() );
  std::string line;

  // create regular expression to get module and per event time
  // prototype: EventCounter_hp_TOP: per event time (ms):    0.414828
  TPRegexp regexp( "(\\w+)_TOP: per event time \\(ms\\):\\s+((?:\\d|.)+)" );

  while( std::getline( in, line ) )
  {
    auto matches = regexp.MatchS( line.c_str() );
    if( matches->GetLast() < 2 ) continue;

    // get module name
    const std::string module( static_cast<TObjString*>(matches->At(1))->GetString().Data() );

    // get time
    const auto time( static_cast<TObjString*>(matches->At(2))->GetString().Atof() );

    // printout
    std::cout << "ParseFile - module: " << module << " time: " << time << std::endl;

    // insert in map
    module_times[module].push_back( time );
  }

}

//________________________________________________________________________________
TCanvas* DrawModuleTime( const std::string& module, const TimeArray& time_array )
{
  if( time_array.empty() ) return nullptr;
  auto cvname = Form( "cv_%s", module.c_str() );
  auto cv = new TCanvas( cvname, cvname, 800, 800 );

  // get min and max values
  auto min = *std::min_element( time_array.begin(), time_array.end() );
  auto max = *std::max_element( time_array.begin(), time_array.end() );

  auto hname = Form( "h_%s", module.c_str() );
  auto h = new TH1F( hname, hname, 100, 0.5*min, 1.5*max );
  h->GetXaxis()->SetTitle( Form( "%s time (ms)", module.c_str() ) );

  for( const auto& time:time_array )
  { h->Fill( time ); }

  h->SetFillStyle(1001);
  h->SetFillColor(kYellow);
  h->Draw();

  Draw::PutText( 0.5, 0.6, Form( "mean: %.1f #pm %.1f (ms)", h->GetMean(), h->GetMeanError() ) );

  return cv;
}

//________________________________________________________________________________
double GetAverage( const TimeArray& time_array )
{ return std::accumulate( time_array.begin(), time_array.end(), (double) 0, [](const double& sum, const double& current ) { return sum+current; } )/time_array.size(); }

//________________________________________________________________________________
TCanvas* DrawPie( const ModuleMap& module_times)
{
  std::vector<double> timelist;
  std::vector<const char*> namelist;
  double sum = 0;
  for( const auto& [name, time_array]:module_times )
  {
    auto average = GetAverage( time_array );
    sum += average;
    namelist.push_back( Form( "%s (%.0fms)", name.c_str(), average));
    timelist.push_back( average );
  }

  // remover labels for too small contributions
  for( int i = 0; i < timelist.size(); ++i )
  { if( timelist[i]/sum < 0.02 ) namelist[i] = ""; }

  auto cv = new TCanvas( "cvpie", "cvpie", 800, 800 );
  auto pie = new TPie("pie", "module time", timelist.size(), &timelist[0] );
  pie->SetLabels( &namelist[0] );
  pie->SetRadius(.2);
  pie->SetTextSize(0.02);
  pie->Draw("<");

  Draw::PutText( 0.3, 0.1, Form("Total: %.1f ms/event", sum  ));

  return cv;

}

//________________________________________________________________________________
void ParseTime()
{
  const TString tag = "_CombinedDataReconstruction_kfp_test3";
  // const TString tag = "_CombinedDataReconstruction_kfp";
  const TString inputFile = Form( "LOG/output%s_*.txt", tag.Data());
  TString pdfFile( Form( "Figures/Time_reco%s.pdf", tag.Data() ) );

  PdfDocument pdfDocument( pdfFile );

  ModuleMap module_times;

  FileManager fileManager( inputFile );
  std::cout << "ParseTime - loaded " << fileManager.GetNFiles() << " files" << std::endl;
  for( const auto& file:fileManager.GetFiles() )
  { ParseFile( file.Data(), module_times ); }

  double sum = 0;
  for( const auto& [name, time_array]:module_times )
  {
    sum += GetAverage( time_array );
    auto cv = DrawModuleTime( name, time_array );
    if( cv ) { pdfDocument.Add( cv ); }
  }

  pdfDocument.Add( DrawPie( module_times ) );

  std::cout << "Time per event: " << sum << std::endl;
}
