#include "AtMAGNEXMap.h"

#include <FairLogger.h>

#include <Math/Point2D.h>
#include <TH2Poly.h>
#include <TString.h>

#include <boost/multi_array/base.hpp>
#include <boost/multi_array/extent_gen.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include <boost/multi_array/subarray.hpp>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>

using XYPoint = ROOT::Math::XYPoint;

AtMAGNEXMap::AtMAGNEXMap(TString pathToPadDefinitionFile) : AtMap()
{
   SetPadParameters(pathToPadDefinitionFile);
   AtPadCoord.resize(boost::extents[fNumberPads][4][2]);
   std::cout << " MAGNEX Map initialized. " << std::endl;
   std::fill(AtPadCoord.data(), AtPadCoord.data() + AtPadCoord.num_elements(), 0);
   std::cout << " MAGNEX Pad Coordinates container initialized. " << std::endl;
}

AtMAGNEXMap::~AtMAGNEXMap() = default;

void AtMAGNEXMap::Dump() {}

void AtMAGNEXMap::GeneratePadPlane()
{
   if (fPadPlane) {
      LOG(error) << " AtMAGNEXMap::GeneratePadPlane Error : The pad plane is already parsed. Skipping generation of pad plane.";
      return;
   }

   for (Int_t iCol = 0; iCol < fColNum; ++iCol) {
      for (Int_t iRow = 0; iRow < fRowNum; ++iRow) {
         AtPadCoord[PadID(iCol, iRow)][0][0] = fColWidth * iCol;
         AtPadCoord[PadID(iCol, iRow)][0][1] = fRowStart + (fRowWidth + fRowSeparation) * iRow;
         AtPadCoord[PadID(iCol, iRow)][1][0] = fColWidth * (iCol + 1);
         AtPadCoord[PadID(iCol, iRow)][1][1] = fRowStart + (fRowWidth + fRowSeparation) * iRow;
         AtPadCoord[PadID(iCol, iRow)][2][0] = fColWidth * (iCol + 1);
         AtPadCoord[PadID(iCol, iRow)][2][1] = fRowStart + (fRowWidth + fRowSeparation) * iRow + fRowWidth;
         AtPadCoord[PadID(iCol, iRow)][3][0] = fColWidth * iCol;
         AtPadCoord[PadID(iCol, iRow)][3][1] = fRowStart + (fRowWidth + fRowSeparation) * iRow + fRowWidth;
      }
   }

   fPadPlane = new TH2Poly();
   for (Int_t iPad = 0; iPad < fNumberPads; ++iPad) {
      Double_t px[] = {AtPadCoord[iPad][0][0], AtPadCoord[iPad][1][0], AtPadCoord[iPad][2][0], AtPadCoord[iPad][3][0],
                       AtPadCoord[iPad][0][0]};
      Double_t py[] = {AtPadCoord[iPad][0][1], AtPadCoord[iPad][1][1], AtPadCoord[iPad][2][1], AtPadCoord[iPad][3][1],
                       AtPadCoord[iPad][0][1]};
      fPadPlane->AddBin(5, px, py);
   }

   fPadPlane->SetName("MAGNEX_Plane");
   fPadPlane->SetTitle("MAGNEX_Plane");
   fPadPlane->ChangePartition(500, 500);

   kIsParsed = true;
   std::cout << " MAGNEX Pad Plane generated. " << std::endl;
}

XYPoint AtMAGNEXMap::CalcPadCenter(Int_t iPad)
{
   if (!kIsParsed) {
      LOG(error) << " AtMAGNEXMap::CalcPadCenter Error : Pad plane has not been generated or parsed.";
      return {-9999, -9999};
   }

   Double_t x = (AtPadCoord[iPad][0][0] + AtPadCoord[iPad][1][0]) / 2.;
   Double_t y = (AtPadCoord[iPad][0][1] + AtPadCoord[iPad][2][1]) / 2.;
   return {x, y};
}

/*
   The following function sets the parameters that define the pad plane. These parameters are the number of
   columns, the number of rows, the width of every column and row, separation of different rows and the
   position of the first ow with respect to the border of the detector.
   Inputs: + TString pathToPadDefinitionFile: the path of the file that contains the values of the parameters
                                              defining the pad plane. The structure of the file must be the
                                              following:
                                              Line 1) #Columns ColumnWidth
                                              Line 2) #Rows    RowWidth     RowSeparation RowStart
   Returns: nothing.
   Status: Finished (I think).
*/
void AtMAGNEXMap::SetPadParameters(TString pathToPadDefinitionFile)
{
   try {
   std::ifstream file;
   file.open(pathToPadDefinitionFile.Data());

   std::cout << "Opened " << pathToPadDefinitionFile.Data() << std::endl;

   std::string line1, line2;
   std::getline(file, line1);
   std::getline(file, line2);

   std::istringstream isxs1(line1);
   std::istringstream isxs2(line2);

   isxs1 >> fColNum >> fColWidth;
   isxs2 >> fRowNum >> fRowWidth >> fRowSeparation >> fRowStart;

   fNumberPads = fColNum * fRowNum;

   std::cout << "Number of columns: " << fColNum << std::endl;
   std::cout << "Number of rows: " << fRowNum << std::endl;

   std::cout << "Column width: " << fColWidth << " mm" << std::endl;

   std::cout << "Row width: " << fRowWidth << " mm" << std::endl;
   std::cout << "Rows separation: " << fRowSeparation << " mm" << std::endl;
   std::cout << "First row position: " << fRowStart << " mm" << std::endl;
   } catch (...) {
      LOG(error) << " AtMAGNEXMap::SetPadParameters Error : " << pathToPadDefinitionFile.Data() << " does not exist, could not be opened or some other error appeared while reading it.";
      return;
   }
}

/*
   The following function simply maps the column and row to the corresponding PadID.
   Inputs: + Int_t StripID: Index corresponding to the strip.
           + Int_t MTHGEMID: Index corresponding to the MTHGEM row.
   Returns: Int_t corresponding with the PadID.
   Status: Finished (I think).
*/
Int_t AtMAGNEXMap::PadID(Int_t iCol, Int_t iRow) {return iCol + iRow * fColNum;}


ClassImp(AtMAGNEXMap)
