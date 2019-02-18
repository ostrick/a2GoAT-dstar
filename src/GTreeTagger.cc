#include "GTreeTagger.h"
#include "GTreeManager.h"

GTreeTagger::GTreeTagger(GTreeManager *Manager)    :
    GTree(Manager, TString("tagger")),
    nTagged(0),
    nDouble(0),
    nChain(0),
    hasEnergy(0),
    hasCalibration(0)
{
    for(Int_t i=0; i<GTreeTagger_MAX; i++)
    {
        taggedChannel[i] = 0;
        taggedTime[i]    = 0;
        taggedEnergy[i]  = 0;
        taggedDouble[i]  = 0;
        doubleChannel[i] = 0;
        doubleRandom[i]  = 0;
        doubleTime[i]    = 0;
        doubleEnergy[i]  = 0;
    }
    for(Int_t i=0; i<352; i++) calibration[i] = 0;
    TaggerPairTimeDiff = new TH1D("TaggerPairTimeDiff","Absolute Time Difference between Tagger Pairs",100,0,20);
}

GTreeTagger::~GTreeTagger()
{

}

void    GTreeTagger::SetBranchAdresses()
{
    if(inputTree->GetBranch("nTagged")) inputTree->SetBranchAddress("nTagged", 	   &nTagged);
    if(inputTree->GetBranch("taggedChannel")) inputTree->SetBranchAddress("taggedChannel", taggedChannel);
    if(inputTree->GetBranch("taggedTime")) inputTree->SetBranchAddress("taggedTime",    taggedTime);
    if(inputTree->GetBranch("taggedEnergy"))
    {
        inputTree->SetBranchAddress("taggedEnergy",  taggedEnergy);
        hasEnergy = true;
    }

}

void    GTreeTagger::SetBranches()
{
    outputTree = inputTree->CloneTree(0);
}

void    GTreeTagger::DecodeDoubles(const Double_t timingRes, const Bool_t decodeChain)
{
    // Local variables for search
    Int_t nPairs = 0;
    Double_t timeDiff = 0;

    Int_t taggedOrder[GTreeTagger_MAX] = {0};
    Int_t pairInd1[GTreeTagger_MAX]    = {0};
    Int_t pairInd2[GTreeTagger_MAX]    = {0};
    Double_t pairTime[GTreeTagger_MAX] = {0};
    Double_t pairDiff[GTreeTagger_MAX] = {0};
    Int_t pairOrder[GTreeTagger_MAX]   = {0};

    // Reset global variables for this event
    nDouble = 0;
    nChain = 0;
    for(Int_t i=0; i<nTagged; i++) taggedDouble[i] = false;

    // Sort channel list into ascending order to make the search easier
    TMath::Sort(nTagged, taggedChannel, taggedOrder, false);

    // Loop over new ordered list
    for(Int_t i=0; i<nTagged; i++)
    {
        for(Int_t j=i+1; j<nTagged; j++)
        {
            // Skip if multi-hit of same channel
            if(taggedChannel[taggedOrder[j]] == taggedChannel[taggedOrder[i]]) continue;

            // Break if beyond neighbor
            if(taggedChannel[taggedOrder[j]] > (taggedChannel[taggedOrder[i]] + 1)) break;

            // Check if neighboring channel is within timing resolution
            timeDiff = (TMath::Abs(taggedTime[taggedOrder[j]]-taggedTime[taggedOrder[i]]));
            TaggerPairTimeDiff->Fill(timeDiff);
            if(timeDiff < timingRes)
            {
                pairInd1[nPairs] = taggedOrder[i];
                pairInd2[nPairs] = taggedOrder[j];
                pairTime[nPairs] = ((taggedTime[taggedOrder[i]]+taggedTime[taggedOrder[j]])/2.0);
                pairDiff[nPairs] = timeDiff;
                nPairs++;
            }

            // Ensure we don't run over the lists
            if(nPairs == GTreeTagger_MAX)
            {
                cout << "Reached maximum (" << nPairs << ") number of tagger pairs" << endl << "     Skipping double decoding for this event!" << endl;
                return;
            }
        }
    }

    // Look for neighboring pairs to form chain
    if(decodeChain)
    {
        // Loop over pairs, looking for multiples that represent a chain
        for(Int_t i=0; i<nPairs; i++)
        {
            for(Int_t j=i+1; j<nPairs; j++)
            {
                // Skip if multi-hit of same channel
                if(taggedChannel[pairInd1[j]] == taggedChannel[pairInd1[i]]) continue;

                // Break if beyond neighbor
                if(taggedChannel[pairInd1[j]] != taggedChannel[pairInd2[i]]) break;

                // Check if neighboring pair is within timing resolution
                timeDiff = (TMath::Abs(pairTime[j]-pairTime[i]));
                if(timeDiff < timingRes)
                {
                    // Increment number of chains if first pair not already part of a chain
                    if(!taggedDouble[pairInd1[i]] && !taggedDouble[pairInd2[i]]) nChain++;

                    // Denote both pairs as belonging to a chain
                    taggedDouble[pairInd1[i]] = true;
                    taggedDouble[pairInd2[i]] = true;

                    taggedDouble[pairInd1[j]] = true;
                    taggedDouble[pairInd2[j]] = true;
                }
            }

            // Skip if already part of a chain
            if(taggedDouble[pairInd1[i]]) continue;
            if(taggedDouble[pairInd2[i]]) continue;

            // Otherwise call this pair a double
            taggedDouble[pairInd1[i]] = true;
            taggedDouble[pairInd2[i]] = true;

            doubleChannel[nDouble] = taggedChannel[pairInd1[i]];
            doubleRandom[nDouble] = doubleChannel[nDouble]+TMath::Nint(gRandom->Rndm());
            doubleTime[nDouble] = pairTime[i];
            doubleEnergy[nDouble] = ((GetTaggedEnergy(pairInd1[i])+GetTaggedEnergy(pairInd2[i]))/2.0);
            nDouble++;
        }
    }
    // Otherwise choose best pairs
    else
    {
        // Reorder pairs by timing difference
        TMath::Sort(nPairs, pairDiff, pairOrder, false);

        // Loop over pairs, removing from list as we go
        for(Int_t i=0; i<nPairs; i++)
        {
            // Skip if already paired as a double
            if(taggedDouble[pairInd1[pairOrder[i]]]) continue;
            if(taggedDouble[pairInd2[pairOrder[i]]]) continue;

            // Otherwise call this pair a double
            taggedDouble[pairInd1[pairOrder[i]]] = true;
            taggedDouble[pairInd2[pairOrder[i]]] = true;

            doubleChannel[nDouble] = taggedChannel[pairInd1[pairOrder[i]]];
            doubleRandom[nDouble] = doubleChannel[nDouble]+TMath::Nint(gRandom->Rndm());
            doubleTime[nDouble] = pairTime[pairOrder[i]];
            doubleEnergy[nDouble] = ((GetTaggedEnergy(pairInd1[pairOrder[i]])+GetTaggedEnergy(pairInd2[pairOrder[i]]))/2.0);
            nDouble++;
        }
    }
    /*
    if((inputTree->GetReadEntry()) < 100){
        cout << "Event " << manager->GetEventNumber() << " - nDouble = " << nDouble << " - nChain = " << nChain << endl << endl;
        for(Int_t i=0; i<nTagged; i++)
        {
            cout << i << "\t" << taggedChannel[i] << "\t" << taggedTime[i] << "   \t" << taggedOrder[i] << "\t" << taggedChannel[taggedOrder[i]] << "\t" << taggedTime[taggedOrder[i]] << "   \t" << taggedDouble[taggedOrder[i]] << endl;
        }
        cout << endl;
    }
    */
}

void    GTreeTagger::SetCalibration(const Int_t nChan, const Double_t *energy)
{
    for(Int_t i=0; i<nChan; i++) calibration[i] = energy[i];
    hasCalibration = true;
}

TLorentzVector  GTreeTagger::GetVectorProtonTarget(const Int_t index)    const
{
    return TLorentzVector(0, 0, taggedEnergy[index], taggedEnergy[index] + (manager->pdgDB->GetParticle("proton")->Mass()*1000));
}
