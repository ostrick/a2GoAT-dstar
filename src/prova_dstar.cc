#include "prova_dstar.h"
#include <stdio.h>
#include <string.h>

prova_dstar::prova_dstar()
{
    TaggerAccScal_fede      = new TH1D("TaggerAccScal_fede","TaggerAccScal_fede",352,0,352);
    TaggerScalers_uncorr    = new TH1D("TaggerScalers_uncorr","TaggerScalers_uncorr",352,0,352);
    LivetimeHist            = new TH1D("LivetimeHisto", "LivetimeHisto",200,0,100);
    TaggerScalers_corr      = new TH1D("TaggerScalers_corr","TaggerScalers_corr",352,0,352);
    //   Doppia_traccia_residua      = new TH1D("Doppia_traccia_residua","Doppia_traccia_residua",10,0,10);
    //   Doppia_traccia_residua_copl      = new TH1D("Doppia_traccia_residua_copl","Doppia_traccia_residua_copl",10,0,10);
    AllScl_fede      = new TH1D("AllScl_fede","AllScl_fede",8705,0,8705);

    //pi0
    time 	= new GH1("time", 	"time", 	1400, -700, 700);
    TaggerChannel    = new GH1 ("TaggerChannel", "TaggerChannel",  352, 0,   352);
    coplanarity = new GH1 ("coplanarity", "coplanarity",    360, 0, 360);
    coplanarity_cut = new GH1 ("coplanarity_cut", "coplanarity_cut",    360, 0, 360);
    //    coplanarity_charged = new GH1 ("coplanarity_charged", "coplanarity_charged",    360, 0, 360);
    //    coplanarity_neutral = new GH1 ("coplanarity_neutral", "coplanarity_neutral",    360, 0, 360);

    theta   = new GH1 ("theta", "theta",    180,    0,  180);
    theta_charged   = new GH1 ("theta_charged", "theta_charged",    180,    0,  180);
    theta_neutral   = new GH1 ("theta_neutral", "theta_neutral",    180,    0,  180);

    PSA_neutral = new GHistBGSub2("PSA_neutral", "PSA_neutral", 45,   5,   50,  250,   0,    1000);
    PSA_charged = new GHistBGSub2("PSA_charged", "PSA_charged", 45,   5,   50,  250,   0,    1000);

    dE_E_CB = new  GHistBGSub2("dE_E_CB", "dE_E_CB", 200,  0,  400, 50, 0, 10);
    dE_E_TAPS = new  GHistBGSub2("dE_E_TAPS", "dE_E_TAPS", 200,  0,  400, 50, 0, 10);

    px_PSA_cut[0] = 0; px_PSA_cut[1] = 37; px_PSA_cut[2] = 41; px_PSA_cut[3] = 42; px_PSA_cut[4] = 43; px_PSA_cut[5] = 43; px_PSA_cut[6] = 0;
    py_PSA_cut[0] = 0; py_PSA_cut[1] =  0; py_PSA_cut[2] = 40; py_PSA_cut[3] = 100; py_PSA_cut[4] = 200; py_PSA_cut[5] = 1000; py_PSA_cut[6] = 1000;
}
prova_dstar::~prova_dstar()
{
}

Bool_t	prova_dstar::Init()
{
    cout << "Initialising physics analysis..." << endl;
    cout << "--------------------------------------------------" << endl << endl;

    if(!InitBackgroundCuts()) return kFALSE;

    //if(!PPhysics::Init()) return kFALSE;
    if(!InitDecodeDoubles()) return kFALSE;
    if(!InitBackgroundCuts()) return kFALSE;
    if(!InitTargetMass()) return kFALSE;
    if(!InitTaggerChannelCuts()) return kFALSE;
    if(!InitTaggerScalers()) return kFALSE;
    if(!InitDisplayScalers()) return kFALSE;
    //if(!PPhysics::Init()) return kFALSE;
    cout << "--------------------------------------------------" << endl;

    return kTRUE;
}

Bool_t	prova_dstar::Start()
{
    //TaggerAccScal->Reset();
    TaggerAccScal_fede->Reset();
    LivetimeHist->Reset();
    TaggerScalers_corr->Reset();
    TaggerScalers_uncorr->Reset();
    AllScl_fede->Reset();
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();
    TraverseValidEvents();
    eventzero = 0;
    conteggio = 0;
    conteggio_N = 0;
    conteggio_C = 0;
    //  cout << "cvssdere" << endl;
    return kTRUE;
    // //  conteggio_check_1 = 0;
    //  conteggio_check_2 = 0;
    //   conteggio_check_3 = 0;
}

void	prova_dstar::ProcessEvent()
{
    Int_t hel = GetTrigger()->GetHelicity();
    //  const Int_t *tp = GetTrigger()->GetTriggerPattern();
    Int_t nerror = 0;

    nerror = GetTrigger()->GetNErrors();


    if(nerror==0)

    {
        /// fill con indice i

        conteggio ++;
        //  cout << "evento numero: " << conteggio << endl;

        if(GetDecodeDoubles()) GetTagger()->DecodeDoubles();
        if(fabs((GetTracks()->GetPhi(0))-(GetTracks()->GetPhi(1)))>140 && fabs((GetTracks()->GetPhi(0))-(GetTracks()->GetPhi(1)))<220) coplanarity_condition = true;
        else coplanarity_condition = false;

        if(GetTracks()->IsNeutral(0) == GetTracks()->IsCharged(1))
        {
            //  cout << "evento numero: " << conteggio << " ha un carico e un neutro" << endl;
            for (Int_t i = 0; i < GetTracks()->GetNTracks(); i++)
            {
                for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
                {
                    if(i==0) {
                        Fillcoplanarity(j, coplanarity, kTRUE);
                    }
                    // cout << "evento numero: " << conteggio << " ha un carico e un neutro, traccia numero: " << i << endl;
                    FillTime_track(i,j,time);
                    Filltheta_track(i, j, theta, kTRUE);
                    FillTaggerChannel_track(i, j, TaggerChannel);
                    if(coplanarity_condition == true)
                    {
                        if(i==0) {
                            //test
                            Fillcoplanarity(j, coplanarity_cut, kTRUE);
                        }
                        if(GetTracks()->IsNeutral(i))
                        {
                            conteggio_N ++;
                            Filltheta_track(i, j, theta_neutral, kTRUE);
                            if(GetTracks()->HasTAPS(i))
                            {
                                FillPSA_track(i, j, PSA_neutral);
                            }
                        }

                        if(GetTracks()->IsCharged(i))
                        {
                            conteggio_C ++;
                            Filltheta_track(i, j, theta_charged, kTRUE);
                            if(GetTracks()->HasCB(i))
                            {
                                FilldE_E_CB_track(i, j, dE_E_CB);
                            }
                            if(GetTracks()->HasTAPS(i))
                            {
                                FillPSA_track(i, j, PSA_charged);
                                FilldE_E_TAPS_track(i, j, dE_E_TAPS);
                            }
                        }
                    }
                    //   cout << GetTracks()->GetNTracks() << " " << *GetNeutralPions()->GetTrackIndex() << "      check_ track  " << endl;

                    // if(troppo_grosso == 1)FillMass(*GetNeutralPions(),  IM_check);
                    // if(troppo_grosso == 0)FillMissingMass(*GetNeutralPions(),  MM_check);
                    //	cout << "fin qui funziona " << endl;
                    //            for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
                    //            {
                    //                // GetNeutralPions()->GetNSubParticles();

                    //                FillTime(*GetTracks(), i, j, time);
                    //                Filltheta(*GetTracks(), i, j, theta, kTRUE);
                    //                FillTaggerChannel(*GetTracks(), i, j, TaggerChannel);
                }
            }
        }
    }
    //    cout << "numero tracce neutre: " << conteggio_N << " e cariche: " << conteggio_C << endl;
}

void	prova_dstar::ProcessScalerRead()
{
    // Fill Tagger Scalers
    FillScalers(GetTC_scaler_min(),GetTC_scaler_max(), TaggerAccScal_fede);

    //static bool Switch = false;
    static Double_t Livetime;
    //  std::cout << scalers->GetScaler(0) << std::endl;
    if(scalers->GetScaler(0))  //Tagger
    {
        //Switch = true;
        //      std::cout << std::endl;
        Livetime = Double_t(scalers->GetScaler(190))/(scalers->GetScaler(191));
    }

    LivetimeHist->Fill(Livetime*100);
    //Double_t corrfactor = Livetime*100;

    for(int i=0; i<=8705; i++)
    {
        Double_t value = Double_t(scalers->GetScaler(i));

        AllScl_fede->Fill(i, value);
    }

    /// TAGGER CONFIGURAZIONE 2015!!!
    /// TAGGER Section AB 0-31
    for(int i=2864; i<=2895; i++)
    {
        Double_t value = Double_t(scalers->GetScaler(i));
        Double_t newvalue = value*Livetime;
        //cout << "value   " << value << "   newvalue   " << newvalue << "   Livetime   " << Livetime << endl;
        Int_t p = i-2864;
        TaggerScalers_uncorr->Fill(p, value);
        TaggerScalers_corr->Fill(i-2864, newvalue);
    }

    /// TAGGER Section CD 31-63
    for(int i=2064; i<=2095; i++)
    {
        Double_t value = Double_t(scalers->GetScaler(i));
        Double_t newvalue = value*Livetime;
        //cout << "value   " << value << "   newvalue   " << newvalue << "   Livetime   " << Livetime << endl;
        Int_t p = i-2064+32;
        TaggerScalers_uncorr->Fill(p, value);
        TaggerScalers_corr->Fill(i-2064+32, newvalue);
    }

    /// TAGGER Section EF 64-95
    for(int i=2576; i<=2607; i++)
    {
        Double_t value = Double_t(scalers->GetScaler(i));
        Double_t newvalue = value*Livetime;
        //cout << "value   " << value << "   newvalue   " << newvalue << "   Livetime   " << Livetime << endl;
        Int_t p = i-2576+64;
        TaggerScalers_uncorr->Fill(p, value);
        TaggerScalers_corr->Fill(i-2576+64, newvalue);
    }

    /// TAGGER Section GH 96-127
    for(int i=2288; i<=2319; i++)
    {
        Double_t value = Double_t(scalers->GetScaler(i));
        Double_t newvalue = value*Livetime;
        //cout << "value   " << value << "   newvalue   " << newvalue << "   Livetime   " << Livetime << endl;
        Int_t p = i-2288+96;
        TaggerScalers_uncorr->Fill(p, value);
        TaggerScalers_corr->Fill(i-2288+96, newvalue);
    }

    /// TAGGER Section IJ 128-159
    for(int i=2000; i<=2031; i++)
    {
        Double_t value = Double_t(scalers->GetScaler(i));
        Double_t newvalue = value*Livetime;
        //cout << "value   " << value << "   newvalue   " << newvalue << "   Livetime   " << Livetime << endl;
        Int_t p = i-2000+128;
        TaggerScalers_uncorr->Fill(p, value);
        TaggerScalers_corr->Fill(i-2000+128, newvalue);
    }

    /// TAGGER Section KL 160-191
    for(int i=2352; i<=2383; i++)
    {
        Double_t value = Double_t(scalers->GetScaler(i));
        Double_t newvalue = value*Livetime;
        //cout << "value   " << value << "   newvalue   " << newvalue << "   Livetime   " << Livetime << endl;
        Int_t p = i-2352+160;
        TaggerScalers_uncorr->Fill(p, value);
        TaggerScalers_corr->Fill(i-2352+160, newvalue);
    }

    /// TAGGER Section MN 192-223
    for(int i=2640; i<=2671; i++)
    {
        Double_t value = Double_t(scalers->GetScaler(i));
        Double_t newvalue = value*Livetime;
        //cout << "value   " << value << "   newvalue   " << newvalue << "   Livetime   " << Livetime << endl;
        Int_t p = i-2640+192;
        TaggerScalers_uncorr->Fill(p, value);
        TaggerScalers_corr->Fill(i-2640+192, newvalue);
    }

    /// TAGGER Section OP 224-255
    for(int i=2032; i<=2063; i++)
    {
        Double_t value = Double_t(scalers->GetScaler(i));
        Double_t newvalue = value*Livetime;
        //cout << "value   " << value << "   newvalue   " << newvalue << "   Livetime   " << Livetime << endl;
        Int_t p = i-2032+224;
        TaggerScalers_uncorr->Fill(p, value);
        TaggerScalers_corr->Fill(i-2032+224, newvalue);
    }

    /// TAGGER Section QR 256-287
    for(int i=2320; i<=2351; i++)
    {
        Double_t value = Double_t(scalers->GetScaler(i));
        Double_t newvalue = value*Livetime;
        //cout << "value   " << value << "   newvalue   " << newvalue << "   Livetime   " << Livetime << endl;
        Int_t p = i-2320+256;
        TaggerScalers_uncorr->Fill(p, value);
        TaggerScalers_corr->Fill(i-2320+256, newvalue);
    }

    /// TAGGER Section ST 288-319
    for(int i=2608; i<=2639; i++)
    {
        Double_t value = Double_t(scalers->GetScaler(i));
        Double_t newvalue = value*Livetime;
        //cout << "value   " << value << "   newvalue   " << newvalue << "   Livetime   " << Livetime << endl;
        Int_t p = i-2608+288;
        TaggerScalers_uncorr->Fill(p, value);
        TaggerScalers_corr->Fill(i-2608+288, newvalue);
    }

    /// TAGGER Section UV 320-351
    for(int i=2896; i<=2927; i++)
    {
        Double_t value = Double_t(scalers->GetScaler(i));
        Double_t newvalue = value*Livetime;
        //cout << "value   " << value << "   newvalue   " << newvalue << "   Livetime   " << Livetime << endl;
        Int_t p = i-2896+320;
        TaggerScalers_uncorr->Fill(p, value);
        TaggerScalers_corr->Fill(i-2896+320, newvalue);
    }

}

Bool_t	prova_dstar::Write()
{
    return GTreeManager::Write();
    GTreeManager::Write(TaggerAccScal_fede);
    GTreeManager::Write(TaggerScalers_uncorr);
    GTreeManager::Write(TaggerScalers_corr);
    GTreeManager::Write(LivetimeHist);
    GTreeManager::Write(AllScl_fede);

}


