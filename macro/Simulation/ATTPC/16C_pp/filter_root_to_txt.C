#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <iostream>
#include <fstream>

void filter_root_to_txt()
{
   // 1. Nombres de archivos
   std::string inputFileName = "./data/attpcsim_Bfield.root";
   std::string outputFileName = "hits_attpcsim_Bfield.txt";

   std::ofstream hitsFile;
   hitsFile.open(outputFileName);

   TFile *file = new TFile(inputFileName.c_str(), "READ");
   
   if (!file || file->IsZombie()) {
      std::cout << "Error: No se pudo abrir " << inputFileName << std::endl;
      return;
   }

   TTreeReader Reader("cbmsim", file);

   // 2. Conectamos las ramas necesarias
   // Nota: No necesitamos leer EventID porque ya sabemos cuál queremos,
   // pero sí necesitamos TrackID para filtrar.
   TTreeReaderArray<Double_t> arrX(Reader, "AtTpcPoint.fX");
   TTreeReaderArray<Double_t> arrY(Reader, "AtTpcPoint.fY");
   TTreeReaderArray<Double_t> arrZ(Reader, "AtTpcPoint.fZ");
   TTreeReaderArray<Double_t> arrELoss(Reader, "AtTpcPoint.fELoss");
   TTreeReaderArray<Int_t> arrTrackID(Reader, "AtTpcPoint.fTrackID"); 

   // 3. Escribimos la cabecera EXACTA que pediste
   hitsFile << "EventID TrackID X Y Z E_Loss\n";

  int eventCounter = 0;

   // 4. Bucle principal: Reader.Next() avanza evento por evento hasta el final
   while (Reader.Next()) {

      // Imprimir progreso cada 100 eventos para saber que no se ha colgado
      if (eventCounter % 100 == 0) {
         std::cout << "Procesando evento: " << eventCounter << std::endl;
      }

      int nHits = arrX.GetSize();

      // 5. Bucle interno: Hits dentro del evento actual
      for (int i = 0; i < nHits; i++) {
         
         // Guardamos: ID_Evento ID_Track X Y Z E
         hitsFile << eventCounter << " " 
                  << arrTrackID[i] << " "
                  << arrX[i] << " " 
                  << arrY[i] << " " 
                  << arrZ[i] << " " 
                  << arrELoss[i] << "\n";
      }

      eventCounter++;
   }

   hitsFile.close();
   file->Close();
   std::cout << "¡Completado! Total eventos procesados: " << eventCounter << std::endl;
   std::cout << "Archivo guardado en: " << outputFileName << std::endl;
}