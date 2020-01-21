#include "DeltaZ_truth.C"

void DeltaZAll_truth()
{

  const auto tags =
  {
    "_flat_full_notpc_single_nz11k_highpt",
    "_flat_full_notpc_single_nominal_highpt",
    "_flat_full_notpc_single_nz3k_highpt",
    "_flat_full_notpc_single_nz1k_highpt",
    "_flat_full_notpc_single_nz500_highpt"

//     "_flat_full_notpc_nz11k_highpt",
//     "_flat_full_notpc_nominal_highpt",
//     "_flat_full_notpc_nz3k_highpt",
//     "_flat_full_notpc_nz1k_highpt",
//     "_flat_full_notpc_nz500_highpt"

//     "_flat_truth_notpc_single_nz11k",
//     "_flat_truth_notpc_single_nominal",
//     "_flat_truth_notpc_single_nz3k",
//     "_flat_truth_notpc_single_nz1k",
//     "_flat_truth_notpc_single_nz500"

//     "_realistic_truth_notpc_nz11k",
//     "_realistic_truth_notpc_nominal",
//     "_realistic_truth_notpc_nz3k",
//     "_realistic_truth_notpc_nz1k",
//     "_realistic_truth_notpc_nz500"
//     "_realistic_truth_notpc_single_nz11k",
//     "_realistic_truth_notpc_single_nominal",
//     "_realistic_truth_notpc_single_nz3k",
//     "_realistic_truth_notpc_single_nz1k",
//     "_realistic_truth_notpc_single_nz500"
  };

  std::vector<TString> outputFiles;
  for( const auto&tag:tags )
  { outputFiles.push_back( DeltaZ_truth( tag ) ); }

  for( const auto& file:outputFiles )
  { std::cout << file << std::endl; }

}
