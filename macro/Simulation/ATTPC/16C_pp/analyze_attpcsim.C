#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TClonesArray.h"
#include <iostream>
#include <fstream>


void analyze_attpcsim()
{
   // 1. Preparamos archivos
   std::string inputFileName = "./data/attpcsim_Bfield.root";
   std::string outputFileName = "hits_attpcsim_event11_scattered.txt";

   std::ofstream hitsFile;
   hitsFile.open(outputFileName);

   TFile *file = new TFile(inputFileName.c_str(), "READ");
   
   // Check de seguridad
   if (!file || file->IsZombie()) {
      std::cout << "Error: No se pudo abrir " << inputFileName << std::endl;
      return;
   }

   TTreeReader Reader("cbmsim", file);

   // Lectura de ramas
   TTreeReaderArray<Double_t> arrX(Reader, "AtTpcPoint.fX");
   TTreeReaderArray<Double_t> arrY(Reader, "AtTpcPoint.fY");
   TTreeReaderArray<Double_t> arrZ(Reader, "AtTpcPoint.fZ");
   TTreeReaderArray<Double_t> arrELoss(Reader, "AtTpcPoint.fELoss");
   
   // IMPORTANTE: Asegúrate de que esta rama existe. 
   // Si el script falla aquí, cambia "fTrackID" por el nombre correcto en tu TTree::Print()
   TTreeReaderArray<Int_t> arrTrackID(Reader, "AtTpcPoint.fTrackID"); 

   // Cabecera del archivo de texto
   hitsFile << "EventID,TrackID,X,Y,Z,E_Loss\n";

   std::cout << "Comenzando analisis..." << std::endl;

   // --- AQUÍ ESTABA EL ERROR ---
   // Definimos la variable eventID antes de usarla
   int eventID = 0; // Empezamos desde el evento 11

   // 4. Bucle sobre eventos
   while (Reader.Next()) {
      
      // Imprimir progreso cada 100 eventos
      if (eventID % 100 == 0) std::cout << "Procesando evento " << eventID << std::endl;

      int nHits = arrX.GetSize();

      // 5. Bucle sobre hits dentro del evento
      for (int i = 0; i < nHits; i++) {
         
         // Escribimos: EventID TrackID X Y Z Energia
         hitsFile << eventID << " " 
                  << arrTrackID[i] << " "
                  << arrX[i] << " " 
                  << arrY[i] << " " 
                  << arrZ[i] << " " 
                  << arrELoss[i] << "\n";
      }

      // Incrementamos el contador de evento al finalizar el evento actual
      eventID++;
   }

   hitsFile.close();
   file->Close();
   std::cout << "¡Exito! Archivo generado correctamente: " << outputFileName << std::endl;
}
