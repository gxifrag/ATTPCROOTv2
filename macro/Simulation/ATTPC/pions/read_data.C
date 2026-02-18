void read_data()
{
   // Ajusta la ruta a tu archivo
   TString dir = getenv("VMCWORKDIR");
   TString file = dir + "/macro/Simulation/ATTPC/pions/data/attpcsim_Bfield.root";
   
   TFile *f = new TFile(file, "READ");
   TTree *t = (TTree*)f->Get("cbmsim");
   
   if(!t) { cout << "ERROR: No hay árbol cbmsim" << endl; return; }
   
   TClonesArray *tracks = new TClonesArray("FairMCTrack");
   t->SetBranchAddress("MCTrack", &tracks);
   
   // ATTPC Points (Hits en el detector)
   TClonesArray *points = new TClonesArray("AtTpcPoint");
   t->SetBranchAddress("AtTpcPoint", &points); // Ojo: a veces se llama "AtTpcPoint" o "ATTPCPoint"
   
   cout << "Total Eventos: " << t->GetEntries() << endl;
   
   for(Int_t i=0; i < t->GetEntries(); i++) {
       t->GetEntry(i);
       Int_t nTracks = tracks->GetEntries();
       Int_t nPoints = (points) ? points->GetEntries() : 0;
       
       cout << "Evento " << i << ": Tracks = " << nTracks << ", Points (Hits) = " << nPoints << endl;
       
       if(nPoints > 0) {
           cout << "   ¡BINGO! Hay datos en este evento. El fallo es el Visualizador." << endl;
           return; 
       }
   }
   cout << "RESULTADO: 0 Points encontrados. El fallo es la Simulación (Geometría o Física)." << endl;
}
