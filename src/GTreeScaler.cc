#include "GTreeScaler.h"

#include <TLeaf.h>


GTreeScaler::GTreeScaler(GTreeManager *Manager)    :
    GTree(Manager, TString("scalers"), kTRUE),
    eventNumber(0),
    eventID(0),
    nmr(0),
    nScalers(0)
{
    for(Int_t i=0; i<GTreeScaler_MAX; i++)
        scalers[i] = 0;
}

GTreeScaler::~GTreeScaler()
{

}

void    GTreeScaler::SetBranchAdresses()
{
    if(inputTree->GetBranch("eventNumber")) inputTree->SetBranchAddress("eventNumber", &eventNumber);
    if(inputTree->GetBranch("eventID")) inputTree->SetBranchAddress("eventID", &eventID);
    if(inputTree->GetBranch("nmr")) inputTree->SetBranchAddress("nmr", &nmr);
    if(inputTree->GetBranch("scalers")) nScalers = inputTree->GetLeaf("scalers")->GetLen();
    if(nScalers<=GTreeScaler_MAX && nScalers>0)
        inputTree->SetBranchAddress("scalers", scalers);
}

void    GTreeScaler::SetBranches()
{
    outputTree = inputTree->CloneTree(0);
}


UInt_t  GTreeScaler::GetScalerEntry(const Int_t number)
{
    if(!IsOpenForInput())
    {
        if(!OpenForInput())
        {
            std::cout << "Can not open Scaler tree in input file." << std::endl;
            return 0;
        }
    }

    for(UInt_t i=0; i<GetNEntries(); i++)
    {
        GetEntry(i);
        if(number<eventNumber)
            return i;
    }
    return GetNEntries();
}


void    GTreeScaler::SetNScalers(const Int_t number)
{
    nScalers = number;
    if(nScalers>GTreeScaler_MAX)
    {
        std::cout << "#ERROR# GTreeScaler::SetNScaler(const Int_t num): Can not handle " << number << " Scalers! Set to max (" << GTreeScaler_MAX << ")." << std::endl;
        nScalers = GTreeScaler_MAX;
    }
}

void    GTreeScaler::Print() const
{
    GTree::Print();

    std::cout << "GTreeScaler: eventNumber->" << eventNumber << " eventID->" << eventID << std::endl;
    std::cout << "             nScalers->" << nScalers << std::endl;
    for(Int_t i=0; i<nScalers; i++)
        std::cout << "Scaler " << i << ": " << scalers[i] << std::endl;
}

void    GTreeScaler::CloneValidEntries()
{
    if(!IsOpenForInput())
    {
        if(!OpenForInput())
        {
            std::cout << "Can not open " << GetName() << " in input file." << std::endl;
            return;
        }
    }
    if(!IsOpenForOutput())
    {
        if(!OpenForOutput())
        {
            std::cout << "Can not create " << GetName() << " in output file." << std::endl;
            return;
        }
    }

    for(UInt_t i=1; i<GetNEntries(); i++)
    {
        GetEntry(i);
        Fill();
    }

}
