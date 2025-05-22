#include <RootUtil/Draw.h>
#include <TFile.h>

R__LOAD_LIBRARY(libRootUtilBase.so)


void RawDataMicromegasCluster_QA()
{

  set_style(false);

  const TString inputfile = "Rootfiles/QA/MicromegasClusterQA-00052202-0000.root";
  auto f = TFile::Open( inputfile );

  const std::vector<TString> hnames =
  {
    "h_MicromegasClusterQA_cluster_multiplicity",
    "h_MicromegasClusterQA_cluster_size",
    "h_MicromegasClusterQA_cluster_charge"
  };

  int icv = 0;
  for( const auto& hname:hnames )
  {
    // get histogram
    auto h2 = static_cast<TH2*>( f->Get(hname) );

    // create canvas
    auto cvname = Form("cv%i", icv++ );
    auto cv = new TCanvas(cvname, "", 800, 800 );
    cv->Divide(4,4);

    // loop over detectors
    for( int i=0; i<16;++i )
    {
      const auto detname = h2->GetXaxis()->GetBinLabel( i+1 );
      cv->cd(i+1);

      auto name = Form("%s%i", hname.Data(), i);
      auto h = h2->ProjectionY( name, i+1, i+1 );
      h->SetFillStyle(1001);
      h->SetFillColor(kYellow);
      h->Draw("hist");

      if( icv==2)
      {
        const double fraction = double(h->GetBinContent(1+1))/h->Integral();
        Draw::PutText( 0.2, 0.7, Form("%s - %.3f", detname, fraction ));
      } else {
        Draw::PutText( 0.2, 0.7, detname );
      }

    }
  }

}
