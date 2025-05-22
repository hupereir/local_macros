#include <RootUtil/Draw.h> 
#include <RootUtil/PdfDocument.h>
#include <RootUtil/Utils.h>

#include <TH3.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

static constexpr int isec_rec = 3;
static constexpr double phi_rec = M_PI*isec_rec/6 + M_PI/12;

//_______________________________________________
TString CompareDistortionCorrections()
{
  
  set_style( false );

  
  const TString tag = "_realistic_micromegas_mm-coarse-new";
  TString pdfFile( Form( "Figures/CopmareDistortionCorrections%s_phi%i.pdf", tag.Data(), isec_rec ) );
  PdfDocument pdfDocument( pdfFile );

  static constexpr int nfiles = 2;
  const std::array<TString, nfiles> inputFiles = {{
    Form( "Rootfiles/Distortions_full%s.root", tag.Data() ),
    "distortion_maps/fluct_average-coarse.root"
  }};
  
  // TFiles
  std::array<std::unique_ptr<TFile>, nfiles> tfiles;
  
  // load histograms
  using histogram_array_t = std::array<TH3*, nfiles>;
  histogram_array_t hentries;
  histogram_array_t hDistortionP_rec;
  histogram_array_t hDistortionR_rec;
  histogram_array_t hDistortionZ_rec;

  for( int i = 0; i < nfiles; ++i ) 
  {
    std::cout << "CompareDistortionCorrection - input file: " << inputFiles[i] << std::endl;
    auto f = TFile::Open( inputFiles[i] );
    tfiles[i].reset( f );
    
    hentries[i] = dynamic_cast<TH3*>(f->Get("hentries_rec"));
    hDistortionP_rec[i] = dynamic_cast<TH3*>(f->Get("hDistortionP_rec"));
    hDistortionR_rec[i] = dynamic_cast<TH3*>(f->Get("hDistortionR_rec"));
    hDistortionZ_rec[i] = dynamic_cast<TH3*>(f->Get("hDistortionZ_rec"));
    
    for( const auto&h: { hentries[i], hDistortionP_rec[i], hDistortionR_rec[i], hDistortionZ_rec[i] } )
    { if(h) h->SetName( Form( "%s_%i", h->GetName(), i ) ); }
    
  }
  
  // get relevant phi bin
  int phi_bin = 0;
  for( const auto h:{ hentries[0], hDistortionP_rec[0], hDistortionR_rec[0], hDistortionZ_rec[0] } )
  {
    if( h )
    {
      phi_bin = h->GetXaxis()->FindBin( phi_rec );
      break;
    }
  }
  std::cout << "CompareDistortionCorrections - phi_rec: " << phi_rec << std::endl;
  std::cout << "CompareDistortionCorrections - phi_bin: " << phi_bin << std::endl;

  // get projections of entries histograms if any
  using projection_array_t = std::array<TH1*, nfiles>;
  projection_array_t entries_proj;
  for( int i = 0; i < nfiles; ++i )
  {
    const auto h = hentries[i];
    if( h )
    {
      h->GetXaxis()->SetRange( phi_bin, phi_bin );    
      entries_proj[i] = h->Project3D( "yz" );
    } else entries_proj[i] = nullptr;
  }
  
  // histogram pairs, for comparison
  using histogram_struct_t = std::tuple<TString, histogram_array_t>;
  std::array<histogram_struct_t, 3> histogram_pairs = {{
    { "r#Delta#phi (cm)", hDistortionP_rec },
    { "#Deltar (cm)", hDistortionR_rec },
    { "#Deltaz (cm)", hDistortionZ_rec } 
  }};

 
  for( const auto& [label, histograms]:histogram_pairs )
  {
    
    auto cv( new TCanvas( "cv", "cv1", 1000, 500 ) );
    cv->Divide( 3, 1 );

    projection_array_t proj;
    for( int i =0; i <nfiles; ++i )
    {
      const auto h = histograms[i];
      h->GetXaxis()->SetRange( phi_bin, phi_bin );
      
      proj[i] = h->Project3D( "yz" );
      proj[i]->SetTitle( "" );
      proj[i]->GetZaxis()->SetTitleOffset( 1.6 );
      proj[i]->GetZaxis()->SetTitle( label );

    }
    
//     cv->cd(1);
//     gPad->SetRightMargin( 0.24 );
//     proj[0]->Draw("colz");
//     
//     cv->cd(2);
//     gPad->SetRightMargin( 0.24 );
//     proj[1]->Add( proj[0], -1 );
//     proj[1]->Draw();

    for( int i =0; i <2; ++i )
    {
      cv->cd( i+1 );
      gPad->SetRightMargin( 0.24 );
      proj[i]->DrawClone("colz" );
    }
    
    cv->cd(3);
    gPad->SetRightMargin( 0.24 );
    auto diff = static_cast<TH1*>( proj[1]->Clone( Form( "%s_diff", proj[1]->GetName() ) ) );
    diff->Add( proj[0], -1 );

    // remove all bins from the diff for which there are not enough entries
    for( int ix = 0; ix < diff->GetNbinsX(); ++ix )
      for( int iy = 0; iy < diff->GetNbinsY(); ++iy )
    {
      for( int i = 0; i < nfiles; ++i )
      {
        if( entries_proj[i] && entries_proj[i]->GetBinContent( ix+1, iy+1 ) < 1000 )
        {
          diff->SetBinContent( ix+1, iy+1, 0 );
          diff->SetBinError( ix+1, iy+1, 0 );
          break;
        }
      }
    }
    
    diff->Draw("colz");
    
    pdfDocument.Add(cv);

  }
  
  return pdfFile;

}
