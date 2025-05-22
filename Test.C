#include <TFile.h>
#include <TTree.h>

#include <g4eval/TrackingEvaluator_hp.h>
R__LOAD_LIBRARY(libg4eval.so)

void Test()
{

  TFile* f = TFile::Open("DST/dst_eval_1k_realistic_full_nominal_new2.root");
  TTree* T = static_cast<TTree*>( f->Get( "T" ) );
  T->Draw( "_tracks._clusters._layer", "_tracks._clusters._layer==0" );

}
