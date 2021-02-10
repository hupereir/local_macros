#include "DeltaZ_truth.C"

void DeltaZAll_truth()
{

  const auto tags =
  {
//
//     "_realistic_truth_notpc_nz11k",
//     "_realistic_truth_notpc_nominal",
//     "_realistic_truth_notpc_nz3k",
//     "_realistic_truth_notpc_nz1k",
//     "_realistic_truth_notpc_nz500",
//     "_realistic_truth_notpc_nz300",
//     "_realistic_truth_notpc_nz100",
//     "_realistic_truth_notpc_nz50",
//     "_realistic_truth_notpc_nz10",
//
//     "_realistic_full_notpc_nz11k",
//     "_realistic_full_notpc_nominal",
//     "_realistic_full_notpc_nz3k",
//     "_realistic_full_notpc_nz1k",
//     "_realistic_full_notpc_nz500",
//     "_realistic_full_notpc_nz300",
//     "_realistic_full_notpc_nz100",
//     "_realistic_full_notpc_nz50",
//     "_realistic_full_notpc_nz10"

    "_realistic_truth_notpc_single_nz11k",
    "_realistic_truth_notpc_single_nominal",
    "_realistic_truth_notpc_single_nz3k",
    "_realistic_truth_notpc_single_nz1k",
    "_realistic_truth_notpc_single_nz500",
    "_realistic_truth_notpc_single_nz300",
    "_realistic_truth_notpc_single_nz100",
    "_realistic_truth_notpc_single_nz50",
    "_realistic_truth_notpc_single_nz10"
  };

  std::vector<TString> outputFiles;
  for( const auto&tag:tags )
  { outputFiles.push_back( DeltaZ_truth( tag ) ); }

  for( const auto& file:outputFiles )
  { std::cout << file << std::endl; }

}
