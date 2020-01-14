#include "DeltaRPhi_truth.C"
#include "DeltaZ_truth.C"

#include <TString.h>

#include <iostream>
#include <vector>

void Resolution_truth()
{

  const TString tag = "_5k_full_notpc_noouter";
  std::vector<TString> files;
  files.push_back( DeltaRPhi_truth( tag ) );
  files.push_back( DeltaZ_truth( tag ) );

  for( const auto& file:files )
  { std::cout << "Resolution_truth - " << file << std::endl; }

}
