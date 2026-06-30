void macro(){
double Sn = 1.2181;

TF1 *fBW = new TF1("fBW", [Sn](double *x, double *p){
   double ER = p[0], Gamma0 = p[1];
   if (x[0] <= Sn) return 0.0;
   return Gamma0*0.159154943/((x[0]-ER)*(x[0]-ER)+Gamma0*Gamma0/4);
}, -1, 8.5, 2);

TF1 *fGauss = new TF1("fGauss","TMath::Gaus(x,0,[0],true)",-1,8.5);

TF1Convolution *conv = new TF1Convolution(fBW, fGauss);
conv->SetRange(-1, 8.5);
conv->SetNofPointsFFT(1000);

TF1 *fTest = new TF1("fTest", *conv, -1, 8.5, conv->GetNpar());
fTest->SetParameters(3.5, 0.04, 0.21);  // ER=3.5, Gamma0=0.3 MeV, sigma=0.05 MeV (valores de prueba)
fTest->SetNpx(2000);
// --- NUEVO: Asignar parámetros a las funciones base ---
fBW->SetParameters(3.5, 0.04); // Parámetros [0] y [1] de fBW
fGauss->SetParameter(0, 0.21); // Parámetro [0] de fGauss

    // --- NUEVO: Cambiar colores para distinguirlas visualmente ---
fTest->SetLineColor(kBlack);   // Convolución en negro
fBW->SetLineColor(kBlue);      // Breit-Wigner en azul
fGauss->SetLineColor(kGreen);  // Gaussiana en verde
fBW->Draw();
fGauss->Draw("same");
fTest->Draw("same");
}
