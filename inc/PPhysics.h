#ifndef __PPhysics_h__
#define __PPhysics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GH1.h"
#include "GConfigFile.h"


class	PPhysics : virtual public GTreeManager
{
private:

	Double_t targetmass;
	Double_t Prompt_low;
	Double_t Prompt_high;
	Double_t Random_low1;
	Double_t Random_high1;
	Double_t Random_low2;
	Double_t Random_high2;
	
	Double_t PvR_ratio;
		
	TLorentzVector beam;
	TLorentzVector target;
	TLorentzVector particle;
	TLorentzVector missingp4;
	
	Double_t time;
	Bool_t 	Prompt;
	Bool_t 	Random;

	Int_t TC_cut_min;
	Int_t TC_cut_max;

	Int_t TC_scaler_min;
	Int_t TC_scaler_max;

    Int_t LT_scaler_clock;
    Int_t LT_scaler_inhib;
    Int_t LT_scaler_tginh;

    TObjArray *scalerHists;
    std::vector<Int_t> nScalerSets;
    std::vector<Int_t> scalerChanL;
    std::vector<Int_t> scalerChanH;
    Int_t nScalerHists;

    Bool_t IsDecodeDoubles;
	
protected:


public:
    PPhysics();
    virtual ~PPhysics();

    virtual Bool_t	Init();
	virtual void	Analyse() {;}
	virtual void	Reconstruct();
    virtual void	ProcessScalerRead();
    virtual Bool_t	Write();

	void	FillMissingMass(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void	FillMissingMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);
	void 	FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning = kFALSE);

    Double_t CalcMissingMass(const GTreeParticle &tree, Int_t particle_index, Int_t tagger_index);
    Double_t CalcMissingEnergy(const GTreeParticle &tree, Int_t particle_index, Int_t tagger_index);
    TLorentzVector CalcMissingP4(const GTreeParticle &tree, Int_t particle_index, Int_t tagger_index);

	void 	FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning = kFALSE, Double_t MM_min = -100000, Double_t MM_max = 100000);
	void 	FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning = kFALSE, Double_t MM_min = -100000, Double_t MM_max = 100000);

	void 	FillTime(const GTreeParticle& tree, GH1* gHist);
	void 	FillTime(const GTreeParticle& tree, Int_t particle_index, GH1* gHist);
	void 	FillTimeCut(const GTreeParticle& tree, GH1* gHist);
	void 	FillTimeCut(const GTreeParticle& tree, Int_t particle_index, GH1* gHist);

	void 	FillMass(const GTreeParticle& tree, GH1* gHist);
	void 	FillMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist);
				
	void	SetTarget(Double_t mass) {target = TLorentzVector(0.,0.,0.,mass);}
	TLorentzVector GetTarget() {return target;}

	void 	SetTC_cut(Int_t cut_min, Int_t cut_max) { TC_cut_min = cut_min; TC_cut_max = cut_max; }
	Int_t 	GetTC_cut_min() { return TC_cut_min;}
	Int_t 	GetTC_cut_max() { return TC_cut_max;}

	void 	SetTC_scalers(Int_t sc_min, Int_t sc_max) { TC_scaler_min = sc_min; TC_scaler_max = sc_max; }
	Int_t 	GetTC_scaler_min() { return TC_scaler_min;}
    Int_t 	GetTC_scaler_max() { return TC_scaler_max;}

    void 	SetLT_scalers(Int_t clock, Int_t inhib, Int_t tginh=0) { LT_scaler_clock = clock; LT_scaler_inhib = inhib; LT_scaler_tginh = tginh; }
    Int_t 	GetLT_scaler_clock() { return LT_scaler_clock;}
    Int_t 	GetLT_scaler_inhib() { return LT_scaler_inhib;}
    Int_t 	GetLT_scaler_tginh() { return LT_scaler_tginh;}

    void 	SetDecodeDoubles(Int_t decode) { IsDecodeDoubles = (Bool_t)decode; }
    Bool_t 	GetDecodeDoubles() { return IsDecodeDoubles;}

	// TH1 routines
	void FillMissingMass(const GTreeParticle& tree, TH1* Hprompt, TH1* Hrandom);
	void FillMissingMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom);
	void FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom);

	void FillTime(const GTreeParticle& tree, TH1* Hist);
	void FillTime(const GTreeParticle& tree, Int_t particle_index, TH1* Hist);
	void FillTimeCut(const GTreeParticle& tree, TH1* Hist);
	void FillTimeCut(const GTreeParticle& tree, Int_t particle_index, TH1* Hist);

	void FillMass(const GTreeParticle& tree, TH1* Hist);
	void FillMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hist);

	void FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max);
	void FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max);

    /// START d* FEDE ////
    void    FillTime_track(Int_t particle_index,Int_t tagger_index, GH1* gHist);

    void    Filltheta_track(Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning);
    void    FillPSA_track(Int_t particle_index, Int_t tagger_index, GHistBGSub2* gHist);

    void    FillTaggerChannel_track(Int_t particle_index,Int_t tagger_index, GH1* gHist);

    void    Fillcoplanarity(Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning);


    void    FilldE_E_CB_track(Int_t track_index, GHistBGSub2* gHist);
    void    FilldE_E_CB_track(Int_t track_index, Int_t tagger_index, GHistBGSub2* gHist);

    void    FilldE_E_TAPS_track(Int_t track_index, GHistBGSub2* gHist);
    void    FilldE_E_TAPS_track(Int_t track_index, Int_t tagger_index, GHistBGSub2* gHist);

    /// end d* FEDE ////

	Double_t CalcCoplanarity(const GTreeParticle& tree1, Int_t particle_index1, const GTreeParticle& tree2, Int_t particle_index2);

    void AddScalerHist(const char* name, Int_t lo, Int_t hi);
    void AddScalerHist(const char* name, Int_t scal, const char* label);
    void FillScalers(Int_t low_scaler_number, Int_t high_scaler_number, TH1* hist, Int_t first_bin=1);
    void GoosyTagger(TH1* hist);
    void GoosyVuprom(TH1* hist);
    void GoosyNewFPD(TH1* hist);
    void GoosyNewFPDRecabled(TH1* hist);

    TH1D* GetScalerHist(Int_t index) {return (TH1D*)scalerHists->At(index);}
    TH1D* GetScalerHist(const char* name) {return (TH1D*)scalerHists->At(scalerHists->IndexOf(scalerHists->FindObject(name)));}

	Bool_t InitBackgroundCuts();
	Bool_t InitTargetMass();
    Bool_t InitTaggerCalibration();
	Bool_t InitTaggerChannelCuts();
	Bool_t InitTaggerScalers();
    Bool_t InitLiveTimeScalers();
    Bool_t InitDisplayScalers();
    Bool_t InitDecodeDoubles();

    Bool_t  RejectTagged(Int_t tagger_index);
    Bool_t  RejectDouble(Int_t tagger_index);
};
#endif
