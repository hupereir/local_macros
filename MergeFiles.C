#include <TString.h>

void MergeFiles()
{

  const int n_output_files = 100;
  const int n_events_per_input = 20;
  const int n_files_per_output = 5;

  // generate output file name
  const TString input_path = "DST/CONDOR_Hijing_Micromegas_50kHz/dst_recoeval_merged_truth_notpc/3-newmaps";
  const TString output_path = "DST/CONDOR_Hijing_Micromegas_50kHz/dst_recoeval_merged_truth_notpc/3-newmaps-merged";

  for( int iout = 0; iout < n_output_files; ++iout )
  {
    const int first_event = n_files_per_output*iout*n_events_per_input;
    const int last_event = first_event + n_files_per_output*n_events_per_input;
    const TString output_file = Form( "%s/dst_recoeval_sHijing_0-12fm_merged_%06i_%06i_truth_notpc.root", output_path.Data(), first_event, last_event );
    std::cout << "writing to " << output_file << std::endl;

    TChain* chain = new TChain( "T" );
    for( int iin = 0; iin < n_files_per_output; ++iin )
    {
      TString input_file;
      if( n_events_per_input <= 1 )
      {
        const int event = first_event + iin;
        input_file = Form( "%s/dst_recoeval_sHijing_0-12fm_merged_%06i_truth_notpc.root", input_path.Data(), event );
      } else {
        const int first_event_input = first_event + iin*n_events_per_input;
        const int last_event_input = first_event_input + n_events_per_input;
        input_file = Form( "%s/dst_recoeval_sHijing_0-12fm_merged_%06i_%06i_truth_notpc.root", input_path.Data(), first_event_input, last_event_input );
      }

      if( std::ifstream( input_file.Data() ) )
      {
        std::cout << "adding " << input_file << std::endl;
        chain->Add( input_file );
      }
    }

    // write
    if( chain )
    {
      std::cout << "MergeFiles - entries: " << chain->GetEntries() << std::endl;
      chain->Merge( output_file );
      delete chain;
    }

  }

}
