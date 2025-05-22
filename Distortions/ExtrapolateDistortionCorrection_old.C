#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <tpccalib/TpcSpaceChargeReconstructionHelper.h>

#include <TH3.h>

R__LOAD_LIBRARY(libRootUtilBase.so)
R__LOAD_LIBRARY(libtpccalib.so)

//_______________________________________________
TString ExtrapolateDistortionMap()
{
  set_style( false );

  // input grid
  // const TString tag = "_full_realistic_micromegas_truth-coarse";
  // const TString tag = "_full_realistic_micromegas_mm-coarse-oldgeom";
  // const TString tag = "_full_realistic_micromegas_mm-coarse-newgeom";
  // const TString tag = "_full_realistic_micromegas_mm-coarse-newgeom2";
  // const TString tag = "_full_realistic_micromegas_mm_fullmap-coarse";
  const TString tag = "_full_realistic_micromegas_truth_genfit_nodistortion";
  const auto inputfilename = Form( "Rootfiles/Distortions%s.root", tag.Data() );
  const auto outputfilename =  Form( "Rootfiles/Distortions%s_extrapolated.root", tag.Data() );

  std::cout << "ExtrapolateDistortionMap - inputfilename: " << inputfilename << std::endl;
  std::cout << "ExtrapolateDistortionMap - outputfilename: " << outputfilename << std::endl;

  auto f = TFile::Open( inputfilename );
  if( !f ) return TString();

  // load input histogram
  #if true
  auto hentries = dynamic_cast<TH3*>(f->Get("hentries_rec"));
  auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hDistortionP_rec"));
  auto hDistortionR_rec = dynamic_cast<TH3*>(f->Get("hDistortionR_rec"));
  auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hDistortionZ_rec"));
  #else
  TH3* hentries = nullptr;
  auto hDistortionP_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionP"));
  auto hDistortionR_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionR"));
  auto hDistortionZ_rec = dynamic_cast<TH3*>(f->Get("hIntDistortionZ"));
  #endif

  // output grid
  RootFile fout( outputfilename );

  for( const auto& h:{hentries, hDistortionP_rec, hDistortionZ_rec, hDistortionR_rec } )
  {
    if( !h ) continue;
    TpcSpaceChargeReconstructionHelper::extrapolate_z(h);
    TpcSpaceChargeReconstructionHelper::extrapolate_phi1(h);
    TpcSpaceChargeReconstructionHelper::extrapolate_phi2(h);
    fout.Add( h );
  }

  return outputfilename;
}
