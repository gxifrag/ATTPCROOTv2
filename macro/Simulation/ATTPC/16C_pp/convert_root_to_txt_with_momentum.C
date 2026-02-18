#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <iostream>
#include <fstream>

void convert_root_to_txt_with_momentum()
{
   // 1. Configuración de archivos
   std::string inputFileName = "./data/attpcsim_Bfield.root";
   std::string outputFileName = "hits_attpcsim_all_events_momentum.txt";

   std::ofstream hitsFile;
   hitsFile.open(outputFileName);

   TFile *file = new TFile(inputFileName.c_str(), "READ");
   
   if (!file || file->IsZombie()) {
      std::cout << "Error: No se pudo abrir " << inputFileName << std::endl;
      return;
   }

   TTreeReader Reader("cbmsim", file);

   // 2. Lectura de ramas (Añadimos Px, Py, Pz)
   TTreeReaderArray<Double_t> arrX(Reader, "AtTpcPoint.fX");
   TTreeReaderArray<Double_t> arrY(Reader, "AtTpcPoint.fY");
   TTreeReaderArray<Double_t> arrZ(Reader, "AtTpcPoint.fZ");
   TTreeReaderArray<Double_t> arrELoss(Reader, "AtTpcPoint.fELoss");
   TTreeReaderArray<Int_t> arrTrackID(Reader, "AtTpcPoint.fTrackID"); 

   // --- NUEVAS RAMAS PARA EL MOMENTO ---
   TTreeReaderArray<Double_t> arrPx(Reader, "AtTpcPoint.fPx");
   TTreeReaderArray<Double_t> arrPy(Reader, "AtTpcPoint.fPy");
   TTreeReaderArray<Double_t> arrPz(Reader, "AtTpcPoint.fPz");

   // 3. Cabecera actualizada (9 columnas en total)
   hitsFile << "EventID TrackID X Y Z E_Loss Px Py Pz\n";

   std::cout << "Procesando eventos con informacion de momento..." << std::endl;

   int eventCounter = 0;

   // 4. Bucle sobre todos los eventos
   while (Reader.Next()) {

      if (eventCounter % 100 == 0) std::cout << "Evento: " << eventCounter << std::endl;

      int nHits = arrX.GetSize();

      // 5. Bucle sobre hits
      for (int i = 0; i < nHits; i++) {
         
         // Guardamos: Evento Track X Y Z Energia Px Py Pz
         hitsFile << eventCounter << " " 
                  << arrTrackID[i] << " "
                  << arrX[i] << " " 
                  << arrY[i] << " " 
                  << arrZ[i] << " " 
                  << arrELoss[i] << " "
                  << arrPx[i] << " "    // Nuevo
                  << arrPy[i] << " "    // Nuevo
                  << arrPz[i] << "\n";  // Nuevo
      }

      eventCounter++;
   }

   hitsFile.close();
   file->Close();
   std::cout << "Archivo generado: " << outputFileName << std::endl;
}
