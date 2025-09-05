/*
    MAGNEX Mapping map
    Authors: + Josema josemanuel.lopez.gonzalez@usc.es
*/

#ifndef ATMAGNEXMAP_H
#define ATMAGNEXMAP_H

#include "AtMap.h"

#include <Math/Point2D.h>
#include <TString.h>

#include <unordered_map>

class AtMAGNEXMap : public AtMap {

public:
   AtMAGNEXMap(TString pathToPadDefinitionFile);
   ~AtMAGNEXMap();

   void Dump() override;                                          // Pure virtual
   void GeneratePadPlane() override;                              // Pure virtual
   ROOT::Math::XYPoint CalcPadCenter(Int_t PadRef) override;      // Pure virtual
   Int_t BinToPad(Int_t binval) override {return binval;};        // Pure virtual

   void SetPadParameters(TString pathToPadDefinitionFile);

   Int_t PadID(Int_t iCol, Int_t iRow);
   Int_t GetColNum() {return fColNum;};
   Int_t GetRowNum() {return fRowNum;};

private:
   Int_t fColNum, fRowNum;
   Double_t fColWidth, fRowWidth, fRowSeparation, fRowStart;

   ClassDefOverride(AtMAGNEXMap, 1);
};

#endif
