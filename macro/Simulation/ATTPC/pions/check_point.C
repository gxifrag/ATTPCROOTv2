void check_point()
{
   // 1. Obtener ruta y archivo
   TString dir = getenv("VMCWORKDIR");
   TString geoFile = dir + "/geometry/ATTPC_H300torr.root"; 
   
   cout << "Cargando geometría desde: " << geoFile << endl;

   // 2. Importar geometría
   if (gGeoManager) delete gGeoManager; // Limpiar si existe uno previo
   TGeoManager::Import(geoFile);

   if (!gGeoManager) {
       cout << "ERROR FATAL: No se pudo cargar el GeoManager." << endl;
       return;
   }

   // 3. Imprimir el Top Volume para asegurar que cargó
   TGeoVolume* top = gGeoManager->GetTopVolume();
   if (top) {
       cout << "Volumen Madre (Top): " << top->GetName() << endl;
   } else {
       cout << "ADVERTENCIA: No se encontró Top Volume. Intentando continuar..." << endl;
   }

   // 4. Sondear el eje Z para encontrar el gas
   // Probamos desde el origen hacia adelante (dirección del haz)
   Double_t z_positions[] = {0., 5., 10., 20., 30., 40., 50., 60., 80., 100.}; // en cm
   
   cout << "\n--- BUSCANDO EL GAS (H_300torr) ---" << endl;
   
   for (int i=0; i<10; i++) {
       Double_t z = z_positions[i];
       
       // Convertimos a cm para el output, pero FindNode usa las unidades del archivo (usualmente cm)
       gGeoManager->SetCurrentPoint(0, 0, z);
       TGeoNode *node = gGeoManager->FindNode(0, 0, z);
       
       if (node) {
           TGeoMaterial *mat = node->GetVolume()->GetMaterial();
           TString matName = mat ? mat->GetName() : "Desconocido";
           TString volName = node->GetVolume()->GetName();
           
           cout << "Z = " << z << " cm  --> Vol: " << volName << " \t| Mat: " << matName;
           
           if (matName.Contains("H_300") || matName.Contains("gas") || matName.Contains("hydrogen")) {
               cout << "  <--- ¡AQUI ESTA EL GAS! (Usar este Z)";
           }
           cout << endl;
       } else {
           cout << "Z = " << z << " cm  --> VACÍO (Fuera del mundo)" << endl;
       }
   }
}
