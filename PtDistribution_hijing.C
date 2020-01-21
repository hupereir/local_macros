#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TDatabasePDG.h>
#include <TCanvas.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

static const auto database = new TDatabasePDG();

//____________________________________________
double get_charge( int pdgcode )
{ return database->GetParticle( pdgcode )->Charge()/3; }

//____________________________________________
// require hits in the first VTX layer and the last TPC layer
bool is_charged_particle( int pdgcode )
{ return database->GetParticle( pdgcode ) != 0; }

//____________________________________________
// require hits in the first VTX layer and the last TPC layer
bool mask_is_valid( int64_t mask )
{
  static const auto create_mask = []()
  {
    int64_t mask = 0;
    mask |= 1LL<<54;
    // mask |= 1LL<<54;
    return mask;
  };

  static const int64_t ref_mask = create_mask();
  return ((mask&ref_mask) == ref_mask);
}

//____________________________________________
void PtDistribution_hijing()
{

  const TString inputFile = "DST/dst_eval_hijing_00000_00100-test.root";
  const TString pdfFile( "Figures/PtDistribution_hijing.pdf" );
  PdfDocument pdfDocument( pdfFile );

  FileManager fileManager( inputFile );
  auto tree = fileManager.GetChain( "T" );
  if( !tree ) return;

  const TString var( "_mc_tracks[]._pt" );
  TCut chargeCut( "is_charged_particle( _mc_tracks[]._pid )" );
  TCut maskCut( "mask_is_valid( _mc_tracks[]._mask )" );

  auto h1 = new TH1F( "h", "", 100, 0, 20 );
  Utils::TreeToHisto( tree, h1->GetName(), var, chargeCut && maskCut, false );

  auto cv = new TCanvas( "cv", "cv", 800, 800 );
  h1->GetXaxis()->SetTitle( "#it{p}_{T} (Gev/#it{c})" );
  h1->SetMarkerStyle(20);
  h1->SetMarkerColor(1);
  h1->SetLineColor(1);
  h1->Draw( "E" );
  gPad->SetLogy( true );

  pdfDocument.Add( cv );

}
