
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TStyle.h>

#include <cmath>
#include <cstdint>
#include <iostream>

/// INTT hit structure
class InttHit : public TObject
{
  
  public:
  
  virtual Bool_t IsEqual(const TObject *obj) const
  {
    const auto objcp = dynamic_cast<const InttHit*>(obj);
    return 
      (module == (objcp->module)) &&
      (chip_id== (objcp->chip_id))&&
      (chan_id== (objcp->chan_id))&&
      (adc == (objcp->adc));
  }
  
  virtual Int_t Compare(const TObject* obj) const
  {
    const auto objcp = dynamic_cast<const InttHit*>(obj);
    if( module != objcp->module) return module <  objcp->module ? -1:1;
    if( chip_id != objcp->chip_id) return chip_id < objcp->chip_id ? -1:1;
    if( chan_id != objcp->chan_id) return chan_id < objcp->chan_id ? -1:1;
    return 0;
  }
  
  virtual Bool_t IsSortable() const 
  { return true; }
  
  int pid = 0;
  int adc = 0;
  int ampl = 0;
  int chip_id = 0;
  int module = 0;
  int chan_id = 0;
  int bco = 0;
  Long64_t bco_full = 0;
  int evt = 0;
  int roc = 0;
  int barrel = 0;
  int layer = 0;
  int ladder = 0;
  int arm = 0;
  int full_fphx = 0;
  int full_roc = 0;
  
  //  protected:
  ClassDef(InttHit, 1)
};

/// INTT event structure
class InttEvent : public TObject
{
  public:
  int evtSeq = 0;
  Long64_t bco = 0;
  int fNhits = 0;
  TClonesArray* fhitArray = nullptr;
  int fNclusters = 0;
  ClassDef(InttEvent, 2)
    
};

namespace
{
  
  // check intt hit validity
  bool skip_intt_hit(InttHit* hit)
  {
    bool isBadHit = (
      (hit->bco==0&&hit->chip_id== 0&&hit->chan_id==  0&&hit->adc==0) 
      || (hit->bco==0&&hit->chip_id==21&&hit->chan_id==126&&hit->adc==6) 
      );
    return isBadHit;
  }
  
  // intt basic clustering
  int get_n_clusters( InttEvent* intt_event )
  {
    int nclusters = 0;
    
    // make sure hits are ordered
    intt_event->fhitArray->Sort();
    
    // build list of hits separated by modules and chips
    vector<InttHit*> vHit[14][26];
    const auto N = intt_event->fNhits;
    for(int ihit = 0; ihit < N; ++ihit)
    {
      auto hit = static_cast<InttHit*>( intt_event->fhitArray->UncheckedAt(ihit) );
      if(skip_intt_hit(hit)){continue;}
      if( hit->module >= 14 || hit->chip_id >= 26 ) continue;
      vHit[hit->module][hit->chip_id].push_back(hit);
    }
    
    // clustering
    for(int imodule = 0; imodule < 14; ++imodule)
      for(int ichip = 0; ichip < 26; ++ichip)
    {
      const auto& vhit_local = vHit[imodule][ichip];
      
      std::vector< std::vector<InttHit*> > vlist;
      std::vector<InttHit*> vgroup;
      InttHit* hit_prev = nullptr;
      for( const auto& hit: vhit_local )
      {
        if( hit_prev && (hit->chan_id-hit_prev->chan_id)!=1 )
        {
          vlist.push_back(vgroup);
          vgroup.clear();
        }
        vgroup.push_back(hit);
        hit_prev=hit;
      }
      
      if( !vgroup.empty() ) vlist.push_back(vgroup);
      nclusters += vlist.size();
    }
    return nclusters;
    
  }
  
  // running intt event
  InttEvent* intt_event = nullptr;

}

//_________________________________________________
void INTT_clusterizer(const int runnumber = 20445, const int intt_id = 0)
{
  
  const TString intt_inputfilename = Form( "/phenix/u/hpereira/sphenix/work/intt/LUSTRE/beam_intt%i-%08i-0000_event_base.root", intt_id, runnumber );
  const TString intt_outputfilename = Form( "Rootfiles/intt_clusters/beam_intt%i-%08i-0000_clusters.root", intt_id, runnumber );

  std::cout << "INTT_Correlation - intt_inputfilename: " << intt_inputfilename << std::endl;
  std::cout << "INTT_Correlation - intt_outputfilename: " << intt_outputfilename << std::endl;
 
  auto intt_tfile = TFile::Open( intt_inputfilename, "READ" );
  auto intt_tree = static_cast<TTree*>( intt_tfile->Get("tree") );
  intt_tree->SetBranchAddress( "event", &intt_event );

  auto intt_tfile_out = TFile::Open( intt_outputfilename, "RECREATE" );
  auto intt_tree_out = new TTree( "tree", "tree" );
  intt_tree_out->Branch( "event", &intt_event );

  // const auto entries = intt_tree->GetEntries();
  const auto entries = 2e5;
  std::cout << "INTT_clusterizer - entries: " << entries << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    if( !(i%1000) ) std::cout << "INTT_clusterizer - entry: " << i << std::endl;
    intt_tree->GetEntry(i);
    intt_event->fNclusters = get_n_clusters(intt_event);
    intt_tree_out->Fill();
  }
  
  // write output tree  
  intt_tfile_out->cd();
  intt_tree_out->Write();
  intt_tfile_out->Close();

}
