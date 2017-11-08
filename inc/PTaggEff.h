#ifndef __PTaggEff_h__
#define __PTaggEff_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"
#include "TF1.h"
#include "TGraph.h"

class	PTaggEff  : public PPhysics
{
private:
    Int_t nTaggerChannels = 80;

    TH1*	TaggerTime;
    GH1*	TaggerAllHits;
    GH1*	TaggerSingles;
    GH1*	TaggerDoubles;
    TH1*	TaggerAccScal;
    TH1*	LiveTimeScal;
    Bool_t  FreeScalers;
    Bool_t  HasAttenuation;
    TGraph* Attenuation;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual Bool_t    Write();
			
public:
    PTaggEff();
    virtual ~PTaggEff();
    virtual Bool_t  Init();
    Bool_t InitFreeScalers();
    Bool_t InitAttenuation();
    TGraph GetAttenuation(TString, Double_t, Double_t);

};
#endif
