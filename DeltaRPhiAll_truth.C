#include "DeltaRPhi_truth.C"

void DeltaRPhiAll_truth()
{

  const auto tags =
  {
//     "_realistic_full_notpc_nphi20k",
//     "_realistic_full_notpc_nominal",
//     "_realistic_full_notpc_nphi5k",
//     "_realistic_full_notpc_nphi2k",
//     "_realistic_full_notpc_nphi1k",
//     "_realistic_full_notpc_nphi500",
//     "_realistic_full_notpc_nphi200",
//     "_realistic_full_notpc_nphi100",
//     "_realistic_full_notpc_nphi50",
//
//     "_realistic_truth_notpc_nphi20k",
//     "_realistic_truth_notpc_nominal",
//     "_realistic_truth_notpc_nphi5k",
//     "_realistic_truth_notpc_nphi2k",
//     "_realistic_truth_notpc_nphi1k",
//     "_realistic_truth_notpc_nphi500",
//     "_realistic_truth_notpc_nphi200",
//     "_realistic_truth_notpc_nphi100",
//     "_realistic_truth_notpc_nphi50"

//     "_realistic_truth_notpc_single_nphi20k",
//     "_realistic_truth_notpc_single_nominal",
    "_realistic_truth_notpc_single_nphi5k",
    "_realistic_truth_notpc_single_nphi2k",
    "_realistic_truth_notpc_single_nphi1k",
    "_realistic_truth_notpc_single_nphi500",
    "_realistic_truth_notpc_single_nphi200",
    "_realistic_truth_notpc_single_nphi100",
    "_realistic_truth_notpc_single_nphi50"

  };

  std::vector<TString> outputFiles;
  for( const auto&tag:tags )
  { outputFiles.push_back( DeltaRPhi_truth( tag ) ); }

  for( const auto& file:outputFiles )
  { std::cout << file << std::endl; }

}
