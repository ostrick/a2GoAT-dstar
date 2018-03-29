#include "PTaggEff.h"

PTaggEff::PTaggEff()
{
    nTaggerChannels = 0;
    scalerRead = 0;
    FreeScalers = true;
    HasAttenuation = false;
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
    if(!InitAttenuation()) return kFALSE;

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

    //GoosyTagger(TaggerAccScal);
    //GoosyVuprom(TaggerAccScal);
    //GoosyNewFPD(TaggerAccScal);
    GoosyNewFPDRecabled(TaggerAccScal);

    return kTRUE;
}

void	PTaggEff::ProcessEvent()
{
    if(nTaggerChannels == 0)
    {
        nTaggerChannels = GetSetupParameters()->GetNTagger();
        TaggerTime = new TH2D("TaggerTime","Tagger - Time",1400,-700,700,nTaggerChannels,0,nTaggerChannels);
        TaggerPreHits = new TH1D("TaggerPreHits","Tagger - Previous Hits",nTaggerChannels,0,nTaggerChannels);
        TaggerCurHits = new TH1D("TaggerCurHits","Tagger - Current Hits",nTaggerChannels,0,nTaggerChannels);
        TaggerAllHits = new GH1("TaggerAllHits","Tagger - All Hits",nTaggerChannels,0,nTaggerChannels);
        TaggerSingles = new GH1("TaggerSingles","Tagger - Single Hits",nTaggerChannels,0,nTaggerChannels);
        TaggerDoubles = new GH1("TaggerDoubles","Tagger - Double Hits",nTaggerChannels,0,nTaggerChannels);
        TaggerHits = new TH2D("TaggerHits","Tagger Hits",GetScalers()->GetNEntries(),0,GetScalers()->GetNEntries(),nTaggerChannels,0,nTaggerChannels);
        TaggerScalers = new TH2D("TaggerScalers","Tagger Scalers",GetScalers()->GetNEntries(),0,GetScalers()->GetNEntries(),nTaggerChannels,0,nTaggerChannels);
        TaggerHitScal = new TH2D("TaggerHitScal","Tagger Hits/Scalers",GetScalers()->GetNEntries(),0,GetScalers()->GetNEntries(),nTaggerChannels,0,nTaggerChannels);
    }

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
    Int_t fs;
    string config = ReadConfig("Free-Running-Scal");
    if(strcmp(config.c_str(), "nokey") == 0)
    {
        cout << "Assuming scalers are free running! " << endl << endl;
        FreeScalers = true;
    }
    else if(sscanf( config.c_str(), "%d\n", &fs) == 1)
    {
        cout << "Setting free running scalers: " << fs << endl << endl;
        FreeScalers = (Bool_t)fs;
    }
    else
    {
        cout << "Free running scalers not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

Bool_t  PTaggEff::InitAttenuation()
{
    char name[256];
    Double_t density = 0;
    Double_t length = 0;
    string config = ReadConfig("Beam-Attenuation");
    if(strcmp(config.c_str(), "nokey") == 0)
    {
        cout << "Assuming no beam attenuation! " << endl << endl;
        HasAttenuation = false;
    }
    else if(sscanf( config.c_str(), "%s %lf %lf\n", name, &length, &density) >= 1)
    {
        cout << "Setting beam attenuation: " << name << ", " << length << " cm, " << density << " g/cm^3" << endl << endl;
        HasAttenuation = true;

        TString target = name;
        std::vector<TGraph> attenuators;
        attenuators.push_back(GetAttenuation("Air",0.001225,900));

        if(target.EqualTo("lh2",TString::kIgnoreCase))
        {
            if(length==0) length=10;
            attenuators.push_back(GetAttenuation("Kapton",0.0125,1.42));
            attenuators.push_back(GetAttenuation("Hydrogen",length,0.071));
        }
        else if(target.EqualTo("fst",TString::kIgnoreCase))
        {
            attenuators.push_back(GetAttenuation("Aluminum",0.002,2.7));
            attenuators.push_back(GetAttenuation("Titanium",0.004,4.506));
            attenuators.push_back(GetAttenuation("Helium",1,0.145));
            attenuators.push_back(GetAttenuation("Butanol",1.2,0.94));
        }
        else if(!target.EqualTo("Air",TString::kIgnoreCase))
        {
            attenuators.push_back(GetAttenuation(target,length,density));
        }

        Attenuation = new TGraph(attenuators[0].GetN());
        Double_t enr1, enr2, att1, att2;

        printf("\nMaterial\t100\t300\t500\t700\t900\t1100\t1300\t1500\n\n");
        for(UInt_t i=0; i<(attenuators.size()); i++)
        {
            printf("%-10s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",attenuators[i].GetTitle(),attenuators[i].Eval(100),attenuators[i].Eval(300),attenuators[i].Eval(500),attenuators[i].Eval(700),attenuators[i].Eval(900),attenuators[i].Eval(1100),attenuators[i].Eval(1300),attenuators[i].Eval(1500));
            for(Int_t j=0; j<(attenuators[i].GetN()); j++)
            {
                Attenuation->GetPoint(j,enr1,att1);
                attenuators[i].GetPoint(j,enr2,att2);
                if(i==0) Attenuation->SetPoint(j,enr2,att2);
                else if(enr1==enr2) Attenuation->SetPoint(j,enr2,att1*att2);
                else
                {
                    printf("\n\nEnergy values in files do not match!\n");
                    return kFALSE;
                }
            }
        }
        printf("\n%-10s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n\n","Total",Attenuation->Eval(100),Attenuation->Eval(300),Attenuation->Eval(500),Attenuation->Eval(700),Attenuation->Eval(900),Attenuation->Eval(1100),Attenuation->Eval(1300),Attenuation->Eval(1500));
    }
    else
    {
        cout << "Beam attenuation not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

TGraph  PTaggEff::GetAttenuation(TString target, Double_t length, Double_t density)
{
    TTree tAtt("tAtt","Beam Attenuation");
    tAtt.ReadFile("configfiles/attenuation/"+target+".txt","E:A");
    Int_t nPnts = tAtt.Draw(Form("E:TMath::Exp(-A*%.6f*%.6f)>>attTemp",length,density),"","goff");
    Double_t *enr = tAtt.GetV1();
    Double_t *att = tAtt.GetV2();
    TGraph gAtt(nPnts,enr,att);
    gAtt.SetTitle(target);
    gAtt.SetLineWidth(2);
    //TH1F *attTemp = (TH1F*)gDirectory->Get("attTemp");
    //delete attTemp;
    gDirectory->Delete("attTemp");
    return gAtt;
}

void	PTaggEff::ProcessScalerRead()
{
    TaggerTime->ProjectionY("TaggerCurHits");
    TaggerCurHits->Add(TaggerPreHits,-1);

    TH1D *TaggerPreScal = (TH1D*)TaggerAccScal->Clone("TaggerPreScal");
    PPhysics::ProcessScalerRead();
    scalerRead++;
    TH1D *TaggerCurScal = (TH1D*)TaggerAccScal->Clone("TaggerCurScal");
    TaggerCurScal->Add(TaggerPreScal,-1);

    GoosyNewFPDRecabled(TaggerCurScal);

    for (Int_t i = 1; i<=nTaggerChannels; i++)
    {
        TaggerHits->SetBinContent(scalerRead,i,TaggerCurHits->GetBinContent(i));
        TaggerScalers->SetBinContent(scalerRead,i,TaggerCurScal->GetBinContent(i));
        TaggerHitScal->SetBinContent(scalerRead,i,(TaggerCurHits->GetBinContent(i)/TaggerCurScal->GetBinContent(i)));
    }

    TaggerTime->ProjectionY("TaggerPreHits");

    delete TaggerPreScal;
    delete TaggerCurScal;
}

Bool_t	PTaggEff::Write()
{    
    TH1D *ScalerAllHits = (TH1D*)TaggerAccScal->Clone("ScalerAllHits");
    TH1D *ScalerSingles = (TH1D*)TaggerAccScal->Clone("ScalerSingles");
    TH1D *ScalerDoubles = (TH1D*)TaggerAccScal->Clone("ScalerDoubles");

    TH1D *ScalerOverlap = new TH1D("ScalerOverlap","Scalers in Overlaps",351,0,351);
    TF1 brem("brem","1/x",0.1,1600);
    Double_t over_int, chan_int;

    TH1D *BmAttenuation = new TH1D("BmAttenuation","Beam Attenuation",nTaggerChannels,0,nTaggerChannels);

    for (Int_t i = 0; i < (GetSetupParameters()->GetNTagger()); i++)
    {
        if(HasAttenuation)
        {
            BmAttenuation->SetBinContent(i+1,Attenuation->Eval(GetSetupParameters()->GetTaggerPhotonEnergy(i+1)));
            BmAttenuation->SetBinError(i+1,0.0001);
        }
        if(i==0) continue;

        if(((TaggerAccScal->GetBinContent(i))==0) || ((TaggerAccScal->GetBinContent(i+1))==0)) continue;

        over_int = brem.Integral(((GetSetupParameters()->GetTaggerPhotonEnergy(i-1))-(0.5*(GetSetupParameters()->GetTaggerEnergyWidth(i-1)))),
                                 ((GetSetupParameters()->GetTaggerPhotonEnergy(i))+(0.5*(GetSetupParameters()->GetTaggerEnergyWidth(i)))));

        if((TaggerAccScal->GetBinContent(i)) < (TaggerAccScal->GetBinContent(i+1)))
        {
            chan_int = brem.Integral(((GetSetupParameters()->GetTaggerPhotonEnergy(i-1))-(0.5*(GetSetupParameters()->GetTaggerEnergyWidth(i-1)))),
                                     ((GetSetupParameters()->GetTaggerPhotonEnergy(i-1))+(0.5*(GetSetupParameters()->GetTaggerEnergyWidth(i-1)))));
            ScalerOverlap->SetBinContent(i,(TaggerAccScal->GetBinContent(i))*over_int/chan_int);
        }
        else
        {
            chan_int = brem.Integral(((GetSetupParameters()->GetTaggerPhotonEnergy(i))-(0.5*(GetSetupParameters()->GetTaggerEnergyWidth(i)))),
                                     ((GetSetupParameters()->GetTaggerPhotonEnergy(i))+(0.5*(GetSetupParameters()->GetTaggerEnergyWidth(i)))));
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

    TaggerHits->Write();
    TaggerScalers->Write();
    TaggerHitScal->Write();
    TaggerTime->Write();

    TH1D *TempAllHits = (TH1D*)TaggerAllHits->GetSum()->GetResult()->GetBuffer()->Clone("TempAllHits");
    TH1D *TempSingles = (TH1D*)TaggerSingles->GetSum()->GetResult()->GetBuffer()->Clone("TempSingles");
    TH1D *TempDoubles = (TH1D*)TaggerDoubles->GetSum()->GetResult()->GetBuffer()->Clone("TempDoubles");

    TH1D *TaggEffAllHits = new TH1D("TaggEffAllHits","Tagging Efficiency - All Hits",nTaggerChannels,0,nTaggerChannels);
    TH1D *TaggEffSingles = new TH1D("TaggEffSingles","Tagging Efficiency - Single Hits",nTaggerChannels,0,nTaggerChannels);
    TH1D *TaggEffDoubles = new TH1D("TaggEffDoubles","Tagging Efficiency - Double Hits",nTaggerChannels,0,nTaggerChannels);

    TaggEffAllHits->Sumw2();
    TaggEffAllHits->Divide(TempAllHits,ScalerAllHits,1,LiveTime);

    TaggEffSingles->Sumw2();
    TaggEffSingles->Divide(TempSingles,ScalerSingles,1,LiveTime);

    TaggEffDoubles->Sumw2();
    TaggEffDoubles->Divide(TempDoubles,ScalerDoubles,1,LiveTime);

    if(HasAttenuation)
    {
        TaggEffAllHits->Divide(BmAttenuation);
        TaggEffSingles->Divide(BmAttenuation);
        TaggEffDoubles->Divide(BmAttenuation);
        BmAttenuation->Write();
        delete BmAttenuation;
    }

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
