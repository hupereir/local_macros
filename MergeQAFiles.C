#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

TString MergeQAFiles( TString tag = TString() )
{
//   if( tag.IsNull() ) tag = "_full_high_occupancy" ;
//   const TString inputFile = Form( "QA/CONDOR%s/qa%s_*.root", tag.Data(), tag.Data() );
//   const TString outputFile = Form( "QA/CONDOR%s/qa%s-merged.root", tag.Data(), tag.Data() );

//   if( tag.IsNull() ) tag = "_full_high_occupancy_notpc" ;
//   const TString inputFile = Form( "QA/CONDOR%s/qa%s_*.root", tag.Data(), tag.Data() );
//   const TString outputFile = Form( "QA/CONDOR%s/qa%s-merged.root", tag.Data(), tag.Data() );
  
//   if( tag.IsNull() ) tag = "_Hijing_Micromegas_100kHz" ;
//   const TString inputFile = Form( "DST/CONDOR%s/qa_output_notpc/qa_output*_notpc.root", tag.Data() );
//   const TString outputFile = Form( "DST/CONDOR%s/qa_output_notpc/qa_output_sHijing_0-12fm_notpc.root", tag.Data() );

  if( tag.IsNull() ) tag = "_Hijing_Micromegas_100kHz" ;
  const TString inputFile = Form( "DST/CONDOR%s/qa_output/qa_output*.root", tag.Data() );
  const TString outputFile = Form( "DST/CONDOR%s/qa_output/qa_output_sHijing_0-12fm.root", tag.Data() );

//   if( tag.IsNull() ) tag = "_Hijing_Micromegas_50kHz" ;
//   const TString inputFile = Form( "DST/CONDOR%s/qa_output_merged/qa_output*.root", tag.Data() );
//   const TString outputFile = Form( "DST/CONDOR%s/qa_output_merged/qa_output_sHijing_0-12fm_merged.root", tag.Data() );

  FileManager( inputFile ).Merge( outputFile );
  return outputFile;
}
