
#include <cstdio>
#include <sstream>
#include <TPRegexp.h>


void ParseRemoteUsage( void )
{

  TPRegexp regexp( "Usr\\s+\\d\\s+(\\d+):(\\d+):(\\d+)" );

  auto h = new TH1F( "time", "time", 500, 3600, 3*3600 );

  const TString tag = "_JPsi_np_pythia8";
  const TString selection = Form( "LOG/log_sphenix%s_*.txt", tag.Data() );
  const TString command = Form( "more %s | grep 'Total Remote Usage'", selection.Data() );
  auto tmp = popen( command.Data(), "r" );
  char line[512];
  while( fgets( line, 512, tmp ) )
  {
    std::cout << line;
    auto matches = regexp.MatchS( line );
    if( matches->GetLast() < 3 ) continue;

    const auto hours =  static_cast<TObjString*>(matches->At(1))->GetString().Atof();
    const auto minutes =  static_cast<TObjString*>(matches->At(2))->GetString().Atof();
    const auto seconds =  static_cast<TObjString*>(matches->At(3))->GetString().Atof();

    const auto time = seconds + minutes*60 + hours * 3600;
    std::cout << "Time: " << time << std::endl;
    h->Fill( time );

  }

  gStyle->SetOptStat(1);
  h->Draw();
  std::cout << "mean: " << h->GetMean() << std::endl;


}
