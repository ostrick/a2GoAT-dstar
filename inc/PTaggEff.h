#ifndef __PTaggEff_h__
#define __PTaggEff_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"

class	PTaggEff  : public PPhysics
{
private:
    TH1*	TaggerTime;
    GH1*	TaggerAllHits;
    GH1*	TaggerSingles;
    GH1*	TaggerDoubles;
    TH1*    ScCorrSingles;
    TH1*    ScCorrDoubles;
    TH1*	TaggerAccScal;
    TH1*	LiveTimeScal;
    Bool_t  FreeScalers;

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

};
#endif
