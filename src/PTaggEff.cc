#include "PTaggEff.h"

PTaggEff::PTaggEff()
{ 
    TaggerTime = new TH2D("TaggerTime","Tagger - Time",1400,-700,700,352,0,352);
    TaggerAllHits = new GH1("TaggerAllHits","Tagger - All Hits",352,0,352);
    TaggerSingles = new GH1("TaggerSingles","Tagger - Single Hits",352,0,352);
    TaggerDoubles = new GH1("TaggerDoubles","Tagger - Double Hits",352,0,352);
    FreeScalers = true;
}

PTaggEff::~PTaggEff()
{
}

Bool_t	PTaggEff::Init()
{
    cout << "Initialising tagging efficiency analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;

    if(!InitBackgroundCuts()) return kFALSE;
    if(!InitFreeScalers()) return kFALSE;

    if(!PPhysics::Init()) return kFALSE;

    TaggerAccScal = GetScalerHist("TaggerAccScal");
    if(!TaggerAccScal)
    {
        cout << "No tagger scaler histogram available" << endl;
        return kFALSE;
    }
    LiveTimeScal = GetScalerHist("LiveTimeScal");
    if(!LiveTimeScal)
    {
        cout << "No live time histogram available" << endl;
        return kFALSE;
    }

    cout << "--------------------------------------------------" << endl;
	return kTRUE;
}

Bool_t	PTaggEff::Start()
{
    if(!IsAcquFile())
    {
        cout << "ERROR: Input File is not an Acqu file." << endl;
        return kFALSE;
    }
    SetAsGoATFile();

    TraverseValidEvents();

    GoosyTagger(TaggerAccScal);
    //GoosyVuprom(TaggerAccScal);

    return kTRUE;
}

void	PTaggEff::ProcessEvent()
{
    if(GetDecodeDoubles()) GetTagger()->DecodeDoubles();

    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
    {
        TaggerTime->Fill(GetTagger()->GetTaggedTime(i),GetTagger()->GetTaggedChannel(i));
        TaggerAllHits->Fill(GetTagger()->GetTaggedChannel(i),GetTagger()->GetTaggedTime(i));
        if(RejectTagged(i)) continue;
        TaggerSingles->Fill(GetTagger()->GetTaggedChannel(i),GetTagger()->GetTaggedTime(i));
        TaggerDoubles->Fill(GetTagger()->GetTaggedChannel(i),GetTagger()->GetTaggedTime(i));
    }
    for (Int_t i = 0; i < GetTagger()->GetNDouble(); i++)
    {
        if(RejectDouble(i)) continue;
        TaggerDoubles->Fill(GetTagger()->GetDoubleRandom(i),GetTagger()->GetDoubleTime(i));
    }
}

Bool_t 	PTaggEff::InitFreeScalers()
{
    Int_t sc1;
    string config = ReadConfig("Free-Running-Scal");
    if(sscanf( config.c_str(), "%d\n", &sc1) == 1)
    {
        cout << "Setting free running scalers: " << sc1 << endl << endl;
        FreeScalers = (Bool_t)sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Free running scalers not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

void	PTaggEff::ProcessScalerRead()
{
    PPhysics::ProcessScalerRead();
}

Bool_t	PTaggEff::Write()
{    
    TH1D *ScalerAllHits = (TH1D*)TaggerAccScal->Clone("ScalerAllHits");
    TH1D *ScalerSingles = (TH1D*)TaggerAccScal->Clone("ScalerSingles");
    TH1D *ScalerDoubles = (TH1D*)TaggerAccScal->Clone("ScalerDoubles");

    TH1D *ScalerOverlap = new TH1D("ScalerOverlap","Scalers in Overlaps",351,0,351);
    TF1 brem("brem","1/x",0.1,1600);
    Double_t over_int, chan_int;

    for (Int_t i = 1; i < (GetSetupParameters()->GetNTagger()); i++)
    {
        if(((TaggerAccScal->GetBinContent(i))==0) || ((TaggerAccScal->GetBinContent(i+1))==0)) continue;

        over_int = brem.Integral(((GetSetupParameters()->GetTaggerPhotonEnergy(i))-(0.5*(GetSetupParameters()->GetTaggerEnergyWidth(i)))),
                                 ((GetSetupParameters()->GetTaggerPhotonEnergy(i+1))+(0.5*(GetSetupParameters()->GetTaggerEnergyWidth(i+1)))));

        if((TaggerAccScal->GetBinContent(i)) < (TaggerAccScal->GetBinContent(i+1)))
        {
            chan_int = brem.Integral(((GetSetupParameters()->GetTaggerPhotonEnergy(i))-(0.5*(GetSetupParameters()->GetTaggerEnergyWidth(i)))),
                                     ((GetSetupParameters()->GetTaggerPhotonEnergy(i))+(0.5*(GetSetupParameters()->GetTaggerEnergyWidth(i)))));
            ScalerOverlap->SetBinContent(i,(TaggerAccScal->GetBinContent(i))*over_int/chan_int);
        }
        else
        {
            chan_int = brem.Integral(((GetSetupParameters()->GetTaggerPhotonEnergy(i+1))-(0.5*(GetSetupParameters()->GetTaggerEnergyWidth(i+1)))),
                                     ((GetSetupParameters()->GetTaggerPhotonEnergy(i+1))+(0.5*(GetSetupParameters()->GetTaggerEnergyWidth(i+1)))));
            ScalerOverlap->SetBinContent(i,(TaggerAccScal->GetBinContent(i+1))*over_int/chan_int);

        }

        ScalerSingles->SetBinContent(i,((ScalerSingles->GetBinContent(i))-(ScalerOverlap->GetBinContent(i))));
        ScalerSingles->SetBinContent(i+1,((ScalerSingles->GetBinContent(i+1))-(ScalerOverlap->GetBinContent(i))));

        ScalerDoubles->SetBinContent(i,((ScalerDoubles->GetBinContent(i))-0.5*(ScalerOverlap->GetBinContent(i))));
        ScalerDoubles->SetBinContent(i+1,((ScalerDoubles->GetBinContent(i+1))-0.5*(ScalerOverlap->GetBinContent(i))));
    }

    Double_t LiveTime = 1;
    if(FreeScalers) LiveTime = ((LiveTimeScal->GetBinContent(2))/(LiveTimeScal->GetBinContent(1)));

    // Write all GH1's and TObjects defined in this class
    if(!(GTreeManager::Write())) return false;

    TH1D *TempAllHits = (TH1D*)TaggerAllHits->GetSum()->GetResult()->GetBuffer()->Clone("TempAllHits");
    TH1D *TempSingles = (TH1D*)TaggerSingles->GetSum()->GetResult()->GetBuffer()->Clone("TempSingles");
    TH1D *TempDoubles = (TH1D*)TaggerDoubles->GetSum()->GetResult()->GetBuffer()->Clone("TempDoubles");

    TH1D *TaggEffAllHits = new TH1D("TaggEffAllHits","Tagging Efficiency - All Hits",352,0,352);
    TH1D *TaggEffSingles = new TH1D("TaggEffSingles","Tagging Efficiency - Single Hits",352,0,352);
    TH1D *TaggEffDoubles = new TH1D("TaggEffDoubles","Tagging Efficiency - Double Hits",352,0,352);

    TaggEffAllHits->Sumw2();
    TaggEffAllHits->Divide(TempAllHits,ScalerAllHits,1,LiveTime);

    TaggEffSingles->Sumw2();
    TaggEffSingles->Divide(TempSingles,ScalerSingles,1,LiveTime);

    TaggEffDoubles->Sumw2();
    TaggEffDoubles->Divide(TempDoubles,ScalerDoubles,1,LiveTime);

    TaggEffAllHits->Write();
    TaggEffSingles->Write();
    TaggEffDoubles->Write();

    ScalerAllHits->Write();
    ScalerSingles->Write();
    ScalerDoubles->Write();
    ScalerOverlap->Write();

    delete ScalerAllHits;
    delete ScalerSingles;
    delete ScalerDoubles;
    delete ScalerOverlap;

    delete TempAllHits;
    delete TempSingles;
    delete TempDoubles;

    delete TaggEffAllHits;
    delete TaggEffSingles;
    delete TaggEffDoubles;

    return true;
}
