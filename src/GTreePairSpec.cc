#include "GTreePairSpec.h"

#include <TLeaf.h>

GTreePairSpec::GTreePairSpec(GTreeManager *Manager)    :
    GTree(Manager, TString("pairSpec"), kTRUE)
{
  for(Int_t i=0; i< 368; i++)
    {
        scalerOpen[i] = 0;
        scalerGated[i] = 0;
	scalerGatedDly[i] = 0;
    }
    
}

GTreePairSpec::~GTreePairSpec()
{

}

void    GTreePairSpec::SetBranchAdresses()
{
    if(inputTree->GetBranch("scalerOpen")) inputTree->SetBranchAddress("scalerOpen", scalerOpen);
    if(inputTree->GetBranch("scalerGated")) inputTree->SetBranchAddress("scalerGated", scalerGated);
    if(inputTree->GetBranch("scalerGatedDly")) inputTree->SetBranchAddress("scalerGatedDly", scalerGatedDly);
}

void    GTreePairSpec::SetBranches()
{
    outputTree = inputTree->CloneTree(0);
}
