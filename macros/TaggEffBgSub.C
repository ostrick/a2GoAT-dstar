void TaggEffBgSub(TString sBeam, TString sBkg1, TString sBkg2="", Bool_t bFreeScal=false){
  TFile fBeam(sBeam,"READ");
  TFile fBkg1(sBkg1,"READ");  
  gROOT->cd();
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);

  TH1D *hBeamAllHits = (TH1D*)fBeam.Get("TaggerAllHits");
  TH1D *hBeamSingles = (TH1D*)fBeam.Get("TaggerSingles");
  TH1D *hBeamDoubles = (TH1D*)fBeam.Get("TaggerDoubles");
  TH1D *hBeamScalAll = (TH1D*)fBeam.Get("ScalerAllHits");
  TH1D *hBeamScalSin = (TH1D*)fBeam.Get("ScalerSingles");
  TH1D *hBeamScalDou = (TH1D*)fBeam.Get("ScalerDoubles");
  TH1D *hBeamLiveTime = (TH1D*)fBeam.Get("LiveTimeScal");

  Double_t dBeamClock = hBeamLiveTime->GetBinContent(hBeamLiveTime->GetXaxis()->FindBin("Clock"));
  Double_t dBeamInhib = hBeamLiveTime->GetBinContent(hBeamLiveTime->GetXaxis()->FindBin("Inhibited"));

  TH1D *hBkg1ScalAll = (TH1D*)fBkg1.Get("ScalerAllHits");
  TH1D *hBkg1ScalSin = (TH1D*)fBkg1.Get("ScalerSingles");
  TH1D *hBkg1ScalDou = (TH1D*)fBkg1.Get("ScalerDoubles");
  TH1D *hBkg1LiveTime = (TH1D*)fBkg1.Get("LiveTimeScal");

  TH1D *hBackScalAll = (TH1D*)hBkg1ScalAll->Clone("hBackScalAll");
  TH1D *hBackScalSin = (TH1D*)hBkg1ScalSin->Clone("hBackScalSin");
  TH1D *hBackScalDou = (TH1D*)hBkg1ScalDou->Clone("hBackScalDou");

  Double_t dBackClock = hBkg1LiveTime->GetBinContent(hBkg1LiveTime->GetXaxis()->FindBin("Clock"));
  Double_t dBackInhib = hBkg1LiveTime->GetBinContent(hBkg1LiveTime->GetXaxis()->FindBin("Inhibited"));
  
  if(sBkg2 != ""){
    TFile fBkg2(sBkg2,"READ");
    gROOT->cd();

    TH1D *hBkg2ScalAll = (TH1D*)fBkg2.Get("ScalerAllHits");
    TH1D *hBkg2ScalSin = (TH1D*)fBkg2.Get("ScalerSingles");
    TH1D *hBkg2ScalDou = (TH1D*)fBkg2.Get("ScalerDoubles");
    TH1D *hBkg2LiveTime = (TH1D*)fBkg2.Get("LiveTimeScal");

    hBackScalAll->Add(hBkg2ScalAll);
    hBackScalSin->Add(hBkg2ScalSin);
    hBackScalDou->Add(hBkg2ScalDou);

    dBackClock += hBkg2LiveTime->GetBinContent(hBkg2LiveTime->GetXaxis()->FindBin("Clock"));
    dBackInhib += hBkg2LiveTime->GetBinContent(hBkg2LiveTime->GetXaxis()->FindBin("Inhibited"));
  }

  TString sOut = sBeam;
  if(sOut.Contains("GoAT")) sOut.ReplaceAll("GoAT","BkgSub");
  else sOut.Prepend("BkgSub_");

  TFile fOut(sOut,"RECREATE");

  TH1D *hEffAllHits = (TH1D*)hBeamAllHits->Clone("TaggEffAllHits");
  TH1D *hEffSingles = (TH1D*)hBeamSingles->Clone("TaggEffSingles");
  TH1D *hEffDoubles = (TH1D*)hBeamDoubles->Clone("TaggEffDoubles");
  TH1D *hEffScalAll = (TH1D*)hBeamScalAll->Clone("TaggEffScalAll");
  TH1D *hEffScalSin = (TH1D*)hBeamScalSin->Clone("TaggEffScalSin");
  TH1D *hEffScalDou = (TH1D*)hBeamScalDou->Clone("TaggEffScalDou");

  hEffScalAll->Sumw2();
  hEffScalSin->Sumw2();
  hEffScalDou->Sumw2();

  if(bFreeScal){
    hEffScalAll->Add(hBackScalAll,(-dBeamClock/dBackClock));
    hEffScalAll->Scale(dBeamInhib/dBeamClock);

    hEffScalSin->Add(hBackScalSin,(-dBeamClock/dBackClock));
    hEffScalSin->Scale(dBeamInhib/dBeamClock);

    hEffScalDou->Add(hBackScalDou,(-dBeamClock/dBackClock));
    hEffScalDou->Scale(dBeamInhib/dBeamClock);
  }
  else{
    hEffScalAll->Add(hBackScalAll,(-dBeamInhib/dBackInhib));
    hEffScalSin->Add(hBackScalSin,(-dBeamInhib/dBackInhib));
    hEffScalDou->Add(hBackScalDou,(-dBeamInhib/dBackInhib));
  }

  hEffAllHits->Divide(hEffScalAll);
  hEffSingles->Divide(hEffScalSin);
  hEffDoubles->Divide(hEffScalDou);

  Double_t dAvg, dMean, dMin, dMax;
  Int_t iMean;
  printf("\n");

  TH1D *hCutAllHits = (TH1D*)hEffAllHits->Clone("TaggEffAllHitsCut");
  TH1D *hProAllHits = new TH1D("ProAllHitsPro","Tagging Efficiency - All Hits",200,0,1);
  TF1 fProAllHits("fProAllHits","gaus",0,352);
  TF1 fAvgAllHits("fAvgAllHits","pol0",0,352);
  TF1 fEffAllHits("fEffAllHits","pol1",0,352);

  dMean = 0;
  iMean = 0;
  for(Int_t i=0; i<hEffAllHits->GetNbinsX(); i++){
    hProAllHits->Fill(hEffAllHits->GetBinContent(i+1));
    if((i > 50) && (i < 250) && (hEffAllHits->GetBinContent(i+1) > 0)){
      dMean += hEffAllHits->GetBinContent(i+1);
      iMean++;
    }
  }
  dAvg = dMean/iMean;
  hProAllHits->Fit("fProAllHits","Q","",(0.5*dAvg),(1.5*dAvg));
  dAvg = fProAllHits.GetParameter(1);
  fAvgAllHits.SetParameter(0,dAvg);

  dMin = (dAvg-5*(fProAllHits.GetParameter(2)));
  dMax = (dAvg+5*(fProAllHits.GetParameter(2)));
  for(Int_t i=0; i<hEffAllHits->GetNbinsX(); i++){
    if((hEffAllHits->GetBinContent(i+1) < dMin) || (hEffAllHits->GetBinContent(i+1) >= dMax)){
      hCutAllHits->SetBinContent(i+1,0);
      hCutAllHits->SetBinError(i+1,0);
    }
  }
  dMax = 0.015*TMath::Nint(100*dAvg);
  hEffAllHits->GetYaxis()->SetRangeUser(0,dMax);
  hCutAllHits->GetYaxis()->SetRangeUser(0,dMax);
  hCutAllHits->Fit("fEffAllHits","Q");
  printf("AllHits: Mean Efficiency = %.3f\tFit Efficiency = %.3f\n",fAvgAllHits.GetParameter(0),fEffAllHits.Eval(176));

  TH1D *hCutSingles = (TH1D*)hEffSingles->Clone("TaggEffSinglesCut");
  TH1D *hProSingles = new TH1D("ProSinglesPro","Tagging Efficiency - Single Hits",200,0,1);
  TF1 fProSingles("fProSingles","gaus",0,352);
  TF1 fAvgSingles("fAvgSingles","pol0",0,352);
  TF1 fEffSingles("fEffSingles","pol1",0,352);

  dMean = 0;
  iMean = 0;
  for(Int_t i=0; i<hEffSingles->GetNbinsX(); i++){
    hProSingles->Fill(hEffSingles->GetBinContent(i+1));
    if((i > 50) && (i < 250) && (hEffSingles->GetBinContent(i+1) > 0)){
      dMean += hEffSingles->GetBinContent(i+1);
      iMean++;
    }
  }
  dAvg = dMean/iMean;
  hProSingles->Fit("fProSingles","Q","",(0.5*dAvg),(1.5*dAvg));
  dAvg = fProSingles.GetParameter(1);
  fAvgSingles.SetParameter(0,dAvg);

  dMin = (dAvg-5*(fProSingles.GetParameter(2)));
  dMax = (dAvg+5*(fProSingles.GetParameter(2)));
  for(Int_t i=0; i<hEffSingles->GetNbinsX(); i++){
    if((hEffSingles->GetBinContent(i+1) < dMin) || (hEffSingles->GetBinContent(i+1) >= dMax)){
      hCutSingles->SetBinContent(i+1,0);
      hCutSingles->SetBinError(i+1,0);
    }
  }
  dMax = 0.015*TMath::Nint(100*dAvg);
  hEffSingles->GetYaxis()->SetRangeUser(0,dMax);
  hCutSingles->GetYaxis()->SetRangeUser(0,dMax);
  hCutSingles->Fit("fEffSingles","Q");
  printf("Singles: Mean Efficiency = %.3f\tFit Efficiency = %.3f\n",fAvgSingles.GetParameter(0),fEffSingles.Eval(176));

  TH1D *hCutDoubles = (TH1D*)hEffDoubles->Clone("TaggEffDoublesCut");
  TH1D *hProDoubles = new TH1D("ProDoublesPro","Tagging Efficiency - Double Hits",200,0,1);
  TF1 fProDoubles("fProDoubles","gaus",0,352);
  TF1 fAvgDoubles("fAvgDoubles","pol0",0,352);
  TF1 fEffDoubles("fEffDoubles","pol1",0,352);

  dMean = 0;
  iMean = 0;
  for(Int_t i=0; i<hEffDoubles->GetNbinsX(); i++){
    hProDoubles->Fill(hEffDoubles->GetBinContent(i+1));
    if((i > 50) && (i < 250) && (hEffDoubles->GetBinContent(i+1) > 0)){
      dMean += hEffDoubles->GetBinContent(i+1);
      iMean++;
    }
  }
  dAvg = dMean/iMean;
  hProDoubles->Fit("fProDoubles","Q","",(0.5*dAvg),(1.5*dAvg));
  dAvg = fProDoubles.GetParameter(1);
  fAvgDoubles.SetParameter(0,dAvg);

  dMin = (dAvg-5*(fProDoubles.GetParameter(2)));
  dMax = (dAvg+5*(fProDoubles.GetParameter(2)));
  for(Int_t i=0; i<hEffDoubles->GetNbinsX(); i++){
    if((hEffDoubles->GetBinContent(i+1) < dMin) || (hEffDoubles->GetBinContent(i+1) >= dMax)){
      hCutDoubles->SetBinContent(i+1,0);
      hCutDoubles->SetBinError(i+1,0);
    }
  }
  dMax = 0.015*TMath::Nint(100*dAvg);
  hEffDoubles->GetYaxis()->SetRangeUser(0,dMax);
  hCutDoubles->GetYaxis()->SetRangeUser(0,dMax);
  hCutDoubles->Fit("fEffDoubles","Q");
  printf("Doubles: Mean Efficiency = %.3f\tFit Efficiency = %.3f\n",fAvgDoubles.GetParameter(0),fEffDoubles.Eval(176));
  printf("\n");

  gROOT->cd();
  c1->Divide(2,3);
  c1->cd(1);
  hEffAllHits->DrawClone();
  fAvgAllHits.DrawClone("same");
  c1->cd(2);
  hCutAllHits->DrawClone();
  c1->cd(3);
  hEffSingles->DrawClone();
  fAvgSingles.DrawClone("same");
  c1->cd(4);
  hCutSingles->DrawClone();
  c1->cd(5);
  hEffDoubles->DrawClone();
  fAvgDoubles.DrawClone("same");
  c1->cd(6);
  hCutDoubles->DrawClone();

  fOut.Write();
  fOut.Close();
}
