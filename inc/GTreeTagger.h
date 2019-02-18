#ifndef __GTreeTagger_h__
#define __GTreeTagger_h__


#include <TLorentzVector.h>
#include "Rtypes.h"
#include "GTree.h"
#include "TRandom3.h"


#define GTreeTagger_MAX 4096


class  GTreeTagger : public GTree
{
private:
    Int_t           nTagged;
    Int_t           taggedChannel[GTreeTagger_MAX];
    Double_t        taggedTime[GTreeTagger_MAX];
    Double_t        taggedEnergy[GTreeTagger_MAX];
    Bool_t          taggedDouble[GTreeTagger_MAX];
    Int_t           nDouble;
    Int_t           nChain;
    Int_t           doubleChannel[GTreeTagger_MAX];
    Int_t           doubleRandom[GTreeTagger_MAX];
    Double_t        doubleTime[GTreeTagger_MAX];
    Double_t        doubleEnergy[GTreeTagger_MAX];
    Bool_t          hasEnergy;
    Bool_t          hasCalibration;
    Double_t        calibration[352];
    TH1*            TaggerPairTimeDiff;

protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTreeTagger(GTreeManager *Manager);
    virtual ~GTreeTagger();

    virtual void    Clear()             {nTagged = 0;}

            Int_t           GetNTagged()                        const	{return nTagged;}
    const	Int_t*          GetTaggedChannel()                  const	{return taggedChannel;}
            Int_t           GetTaggedChannel(const Int_t index) const	{return taggedChannel[index];}
    const	Double_t*       GetTaggedTime()                     const	{return taggedTime;}
            Double_t        GetTaggedTime(const Int_t index)    const	{return taggedTime[index];}
    const	Double_t*       GetTaggedEnergy()                   const	{return taggedEnergy;}
            Double_t        GetTaggedEnergy(const Int_t index)	const	{if(hasEnergy && !hasCalibration) return taggedEnergy[index]; return calibration[taggedChannel[index]];}
    const	Bool_t*         GetTaggedDouble()                   const	{return taggedDouble;}
            Bool_t          GetTaggedDouble(const Int_t index)  const	{return taggedDouble[index];}
            Int_t           GetNDouble()                        const	{return nDouble;}
            Int_t           GetNChain()                         const	{return nChain;}
    const	Int_t*          GetDoubleChannel()                  const	{return doubleChannel;}
            Int_t           GetDoubleChannel(const Int_t index) const	{return doubleChannel[index];}
            Int_t           GetDoubleRandom(const Int_t index)  const	{return doubleRandom[index];}
    const	Double_t*       GetDoubleTime()                     const	{return doubleTime;}
            Double_t        GetDoubleTime(const Int_t index)    const	{return doubleTime[index];}
    const	Double_t*       GetDoubleEnergy()                   const	{return doubleEnergy;}
            Double_t        GetDoubleEnergy(const Int_t index)  const	{return doubleEnergy[index];}
            void            DecodeDoubles(const Double_t timingRes = 3.0, const Bool_t decodeChain = true);
            Bool_t          HasEnergy()                         const   {return hasEnergy;}
            Bool_t          HasCalibration()                    const   {return hasCalibration;}
            void            SetCalibration(const Int_t nChan, const Double_t *energy);
    TLorentzVector          GetVector(const Int_t index)        const   {return TLorentzVector(0, 0, taggedEnergy[index], taggedEnergy[index]);}
    TLorentzVector          GetVectorProtonTarget(const Int_t index)    const;
};

#endif
