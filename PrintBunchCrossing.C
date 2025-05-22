
#include <ffaobjects/EventHeader.h>
#include <ffaobjects/EventHeaderv1.h>
#include <ffaobjects/EventHeaderv2.h>

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

R__LOAD_LIBRARY(libffaobjects.so)
R__LOAD_LIBRARY(libphg4hit.so)
R__LOAD_LIBRARY(libg4testbench.so)
R__LOAD_LIBRARY(libRootUtilBase.so)

void PrintBunchCrossing( int index = 0 )
{
  
  int first_event = 1000*index;
  int last_event = first_event+1000;
  auto filename = Form( "DST/CONDOR_Hijing_Micromegas_50kHz/G4Hits_merged/G4Hits_sHijing_0-12fm_merged_%06i_%06i.root", first_event, last_event );
  std::cout << "PrintBunchCrossing - index: " << index << " filename: " << filename << std::endl;             
  
  auto f = new FileManager( filename );
  auto tree = f->GetChain( "T" );
  
//   auto f = new TFile( "DST/CONDOR_Hijing_Micromegas_50kHz/G4Hits_merged/G4Hits_sHijing_0-12fm_merged_000000_001000.root" );
//   auto tree = static_cast<TTree*>(f->Get("T"));

   std::cout << "PrintBunchCrossing - entries: " << tree->GetEntries() << std::endl;

  // output file 
  std::ofstream out( Form( "timestamps_trigger_50kHz_%06i_%06i.txt", first_event, last_event ) );
  out << "// bunchcrossin id; time (ns)" << std::endl;
  out << "// assuming 106ns between bunches" << std::endl;
  out << "// file: " << filename << std::endl;
  auto header = new EventHeaderv2;
  tree->SetBranchAddress( "DST#EventHeader", &header );
  for( int i = 0; i < tree->GetEntries(); ++i )
  {
    if( !(i%100) ) std::cout << "PrintBunchCrossing - entry: " << i << std::endl; 
    tree->GetEntry( i );
    auto bunchCrossing = header->get_BunchCrossing();
    auto time = bunchCrossing*106;
    out << bunchCrossing << " " << time << std::endl;
  }
  out << "done" << std::endl;
  out.close();
  
}
