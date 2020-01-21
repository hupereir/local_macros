#include "DeltaRPhi_truth.C"

void DeltaRPhiAll_truth()
{

  const auto tags =
  {
    "_flat_full_notpc_single_nphi20k_highpt",
    "_flat_full_notpc_single_nominal_highpt",
    "_flat_full_notpc_single_nphi5k_highpt",
    "_flat_full_notpc_single_nphi2k_highpt",
    "_flat_full_notpc_single_nphi1k_highpt"

//     "_flat_full_notpc_nphi20k_highpt",
//     "_flat_full_notpc_nominal_highpt",
//     "_flat_full_notpc_nphi5k_highpt",
//     "_flat_full_notpc_nphi2k_highpt",
//     "_flat_full_notpc_nphi1k_highpt"
  };

  std::vector<TString> outputFiles;
  for( const auto&tag:tags )
  { outputFiles.push_back( DeltaRPhi_truth( tag ) ); }

  for( const auto& file:outputFiles )
  { std::cout << file << std::endl; }

}
