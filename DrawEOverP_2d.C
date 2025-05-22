#include <RootUtil/PdfDocument.h>

R__LOAD_LIBRARY(libRootUtilBase.so)

void DrawEOverP_2d()
{
  using histogram_pair_t = std::pair<std::string, std::string>;
  using histogram_list_t = std::vector<histogram_pair_t>;

  histogram_list_t histogram_info_list = {
    { "h2d__single_electron","electron" },
    { "h2d__single_positron","positron" },
		{ "h2d__single_piplus","#pi^{+}" },
		{ "h2d__single_piminus","#pi^{-}" },
		{ "h2d__single_kplus","K^{+}" },
		{ "h2d__single_kminus","K^{+}" },
		{ "h2d__single_proton","proton" },
		{ "h2d__single_antiproton","anti-proton" }
  };

  auto f = TFile::Open( "Rootfiles/e_over_p_2d_all.root" );
  PdfDocument pdfDocument( "Figures/e_over_p_2d_all.pdf" );

  std::vector<TH1*> histograms;

  for( const auto& [hname,particle]:histogram_info_list )
  {
    auto h = static_cast<TH1*>(f->Get(hname.c_str()) );
    histograms.push_back(h);
  }

  // sum histograms
  auto h_hadrons = static_cast<TH1*>(histograms[2]->Clone("h_hadrons"));
  h_hadrons->Reset();
  h_hadrons->SetTitle( "hadrons" );
  for( const int& index:{2, 3, 4, 5, 6, 7} )
  {
    h_hadrons->Add(histograms[index]);
  }

  auto h_electrons = static_cast<TH1*>(histograms[0]->Clone("h_electrons"));
  h_electrons->Reset();
  h_electrons->SetTitle( "electrons" );
  h_electrons->Add( histograms[0] );
  h_electrons->Add( histograms[1] );

  // draw
  auto cv = new TCanvas( "cv", "cv", 1200, 600 );
  cv->Divide(2,1);
  cv->cd(1);
  h_hadrons->Draw();

  cv->cd(2);
  h_electrons->Draw();

  pdfDocument.Add(cv);
}
