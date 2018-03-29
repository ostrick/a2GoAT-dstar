#ifndef __GTreeScaler_h__
#define __GTreeScaler_h__


#include "Rtypes.h"
#include "GTree.h"


#define GTreeScaler_MAX 16384

class  GTreeScaler : public GTree
{
private:
    Int_t		eventNumber;
    Int_t		eventID;
    Float_t    nmr;
    UInt_t		scalers[GTreeScaler_MAX];

    Int_t		nScalers;

protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTreeScaler(GTreeManager *Manager);
    virtual ~GTreeScaler();

    virtual void        Clear()     {}
            void        CloneValidEntries();
            Int_t 		GetEventID()        		const	{return eventID;}
            Int_t 		GetEventNumber()        	const	{return eventNumber;}
            Float_t    GetNMR()                    const   {return nmr;}
            Int_t 		GetNScalers()    			const	{return nScalers;}
    const	UInt_t*		GetScaler()                 const	{return	scalers;}
            UInt_t		GetScaler(const Int_t index)const	{return	scalers[index];}
            UInt_t		GetScalerEntry(const Int_t event_number);
            void 		SetNScalers(const Int_t num);
    virtual void        Print() const;
};


#endif
