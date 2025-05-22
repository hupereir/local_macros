R__LOAD_LIBRARY(libmicromegas.so)

void TestNoise(  int runNumber = 57605 )
{

  // const TString inputFile = Form( "Rootfiles/MicromegasCombinedDataEvaluation_%08i-0000.root", runNumber );
  const TString inputFile = "Rootfiles/MicromegasCombinedDataEvaluation_00057605-0000.root";
  auto tfile = TFile::Open( inputFile );
  auto tree = static_cast<TTree*>( tfile->Get( "T" ) );

  tree->Print();

  MicromegasCombinedDataEvaluation::Container* container = nullptr;
  tree->SetBranchAddress( "Event", &container );
  // tree->SetBranchAddress( "DST#EVAL#MicromegasCombinedDataEvaluation::Container", &container );
}
