#include "ATEvent.hh"
#include "ATHit.hh"
#include "ATHoughSpace.hh"
#include "ATHoughSpaceCircle.hh"
#include "ATHoughSpaceLine.hh"
#include "ATPad.hh"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRunAna.h"
#include "MCsrc_cuda.hh"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <ios>
#include <iostream>
#include <istream>
#include <limits>
#include <map>
#include <vector>

Int_t main()
{

   gSystem->Load("libATTPCReco.so");

   FairRunAna *run = new FairRunAna(); // Forcing a dummy run

   TString workdir = getenv("VMCWORKDIR");
   TString FileNameHead = "output";
   TString FilePath = workdir + "/macro/Unpack_GETDecoder2/";
   TString FileNameTail = ".root";
   TString FileName = FilePath + FileNameHead + FileNameTail;

   std::cout << " Opening File : " << FileName.Data() << std::endl;
   TFile *file = new TFile(FileName.Data(), "READ");

   TTree *tree = (TTree *)file->Get("cbmsim");
   Int_t nEvents = tree->GetEntries();
   std::cout << " Number of events : " << nEvents << std::endl;

   TTreeReader Reader1("cbmsim", file);
   TTreeReaderValue<TClonesArray> eventArray(Reader1, "ATEventH");
   TTreeReaderValue<TClonesArray> houghArray(Reader1, "ATHough");

   while (Reader1.Next()) {

      ATEvent *event = (ATEvent *)eventArray->At(0);
      Int_t nHits = event->GetNumHits();
      std::vector<ATHit> *hitArray = event->GetHitArray(); // Not working!
      event->GetHitArrayObj();
      // std::cout<<event->GetHitPadMult(0)<<std::endl;
      // std::cout<<event->GetEventID()<<std::endl;
      hitArray->size();

      std::vector<ATHit *> *hitbuff = new std::vector<ATHit *>; // Working!

      // std::vector<ATEvent*> test;
      // test.push_back(event);

      for (Int_t iHit = 0; iHit < nHits; iHit++) {
         ATHit hit = event->GetHit(iHit);
         TVector3 hitPos = hit.GetPosition();
         hitbuff->push_back(&hit);
      }

      // std::cout<<hitbuff->size()<<std::endl;

      ATHoughSpaceCircle *fHoughSpaceCircle = dynamic_cast<ATHoughSpaceCircle *>(houghArray->At(0));
      // if(!fHoughSpaceCircle) std::cout<<" Warning : Failed casting "<<std::endl;
      std::cout << fHoughSpaceCircle->GetYCenter() << std::endl;
   }

   // #pragma omp parallel for ordered schedule(dynamic,1)
   // for(Int_t i=0;i<100;i++)std::cout<<" Hello ATTPCer! "<<std::endl;
   return 0;
}
