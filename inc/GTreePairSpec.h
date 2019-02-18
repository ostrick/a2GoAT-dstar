#ifndef __GTreePairSpec_h__
#define __GTreePairSpec_h__


#include "Rtypes.h"
#include "GTree.h"

#define GTreePairSpec_MAX 4096

class  GTreePairSpec : public GTree
{
private:
  
    Double_t    scalerOpen[368];
    Double_t	scalerGated[368];
    Double_t	scalerGatedDly[368];

protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTreePairSpec(GTreeManager *Manager);
    virtual ~GTreePairSpec();

    virtual void        Clear()     {}
     const  Double_t*	GetScalerOpen()               const	{return scalerOpen;}
     const  Double_t	GetScalerOpen(Int_t channel)  const	{return scalerOpen[channel];}
     const  Double_t*	GetScalerGated()               const	{return scalerGated;}
     const  Double_t	GetScalerGated(Int_t channel)  const	{return scalerGated[channel];}
     const  Double_t*	GetScalerGatedDly()               const	{return scalerGatedDly;}
     const  Double_t	GetScalerGatedDly(Int_t channel)  const	{return scalerGatedDly[channel];}

};


#endif
