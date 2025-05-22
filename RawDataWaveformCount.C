#include <RootUtil/Draw.h>
#include <RootUtil/FileManager.h>
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>
R__LOAD_LIBRARY(libRootUtilBase.so)

TString RawDataWaveformCount()
{
  // map nunit to file name
  using filename_pair_t = std::pair<int, TString>;
  using filename_map_t = std::list<filename_pair_t>;
//   filename_map_t filename_map = 
//   { 
//     {1, "DST/MicromegasRawDataEvaluation-00013039-0000-full.root" },
//     {16, "DST/MicromegasRawDataEvaluation-00013042-0000-full.root" },
//     {32, "DST/MicromegasRawDataEvaluation-00013043-0000-full.root" },
//     {48, "DST/MicromegasRawDataEvaluation-00013044-0000-full.root" },
//     {64, "DST/MicromegasRawDataEvaluation-00013045-0000-full.root" }
//   };
//   
//   PdfDocument pdfDocument( "Figures/RawDataWaveformCount.pdf" );

//   filename_map_t filename_map = 
//   { 
//     {1, "DST/MicromegasRawDataEvaluation-00014073-0000-full.root" },
//     {8, "DST/MicromegasRawDataEvaluation-00014079-0000-full.root" },
//     {16, "DST/MicromegasRawDataEvaluation-00014074-0000-full.root" },
//     {32, "DST/MicromegasRawDataEvaluation-00014075-0000-full.root" },
//     {48, "DST/MicromegasRawDataEvaluation-00014076-0000-full.root" },
//     {64, "DST/MicromegasRawDataEvaluation-00014077-0000-full.root" }
//   };
//   
//   PdfDocument pdfDocument( "Figures/RawDataWaveformCount-new.pdf" );

  
//   - 14085 NUnit 16 <- after d.reset()
// - 14093 NUnit 1
// - 14094 NUnit 1
// - 14086 NUnit 8
// - 14087 NUnit 16
// - 14088 NUnit 24
// - 14089 NUnit 32
// - 14090 NUnit 40
// - 14091 NUnit 48
// - 14092 NUnit 64
// - 14095 NUnit 16 <- short
// - 14096 NUnit 16
// - 14097 NUnit 8
// - 14098 NUnit 1
// - 14099 NUnit 16

  filename_map_t filename_map = 
  { 
    {1, "DST/MicromegasRawDataEvaluation-00014093-0000-full.root" },
    {1, "DST/MicromegasRawDataEvaluation-00014094-0000-full.root" },
    {1, "DST/MicromegasRawDataEvaluation-00014098-0000-full.root" },
    {8, "DST/MicromegasRawDataEvaluation-00014086-0000-full.root" },
    {8, "DST/MicromegasRawDataEvaluation-00014097-0000-full.root" },
    {16, "DST/MicromegasRawDataEvaluation-00014085-0000-full.root" },
    {16, "DST/MicromegasRawDataEvaluation-00014087-0000-full.root" },
    {16, "DST/MicromegasRawDataEvaluation-00014096-0000-full.root" },
    {16, "DST/MicromegasRawDataEvaluation-00014099-0000-full.root" },
    {24, "DST/MicromegasRawDataEvaluation-00014088-0000-full.root" },
    {32, "DST/MicromegasRawDataEvaluation-00014089-0000-full.root" },
    {40, "DST/MicromegasRawDataEvaluation-00014090-0000-full.root" },
    {48, "DST/MicromegasRawDataEvaluation-00014091-0000-full.root" },
    {64, "DST/MicromegasRawDataEvaluation-00014092-0000-full.root" }
  };
  
  PdfDocument pdfDocument( "Figures/RawDataWaveformCount-new.pdf" );  gStyle->SetOptStat(111111);
  
  std::array<int,2> color = { 1,2 };
  for( const auto& [nunit, filename]:filename_map )
  {
    // open file, get tree
    FileManager f( filename );
    auto tree = f.GetChain( "T" );
    
    // loop over packets
    std::array<TH1*,2> h;
    for( int ipacket = 0; ipacket<2; ++ipacket )
    {
      const auto var = Form("n_waveform[%i]", ipacket);
      const auto hname = Form( "h_%i_%i", nunit, ipacket );
      h[ipacket] = new TH1I( hname, hname, 100, 0, 4096 );
      h[ipacket]->SetLineColor( color[ipacket] );
      h[ipacket]->GetXaxis()->SetTitle( "n_{wf}" );
      Utils::TreeToHisto( tree, h[ipacket]->GetName(), var, TCut(), false );
    
      std::cout << "RawDataWaveformCount - packet: " << ipacket << " overflow: " << h[ipacket]->GetBinContent( h[ipacket]->GetNbinsX()+1 ) << std::endl;

    }
    
    // draw
    const auto cvname = Form( "cv%i", nunit );
    auto cv = new TCanvas( cvname, cvname, 1200, 900 );
    h[0]->SetMaximum( 1.2*std::max(h[0]->GetMaximum(), h[1]->GetMaximum()) );
    h[0]->Draw();
    h[1]->Draw("same");
    gPad->SetLogy(true);
    pdfDocument.Add(cv);
        
  }
  
  return pdfDocument.Filename();

}
