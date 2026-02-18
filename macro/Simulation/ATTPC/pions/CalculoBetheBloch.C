#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <iostream>

// Constantes físicas
const Double_t K = 0.307075;     // MeV cm^2 / mol
const Double_t me = 0.511;       // Masa electrón (MeV)
const Double_t mp = 938.27;      // Masa protón (MeV)
const Double_t rho_STP = 8.375e-5; // Densidad H2 a 1 atm (760 Torr) en g/cm3 (aprox)

// Propiedades del Material (Hidrógeno Gas)
const Double_t Z = 1;            // Número atómico
const Double_t A = 1.0079;       // Masa atómica (g/mol)
const Double_t I = 19.2e-6;      // Potencial de excitación medio (MeV) para H2

// Función para calcular dE/dx (Stopping Power) en MeV/cm
Double_t BetheBloch(Double_t energyMeV, Double_t density) {
    // Kinematics
    Double_t gamma = (energyMeV + mp) / mp;
    Double_t beta = TMath::Sqrt(1 - 1/(gamma*gamma));
    
    if (beta <= 0) return 0;

    Double_t me_c2 = me; 
    
    // Término principal de Bethe-Bloch
    Double_t term1 = (K * Z * Z * Z / A) * (1 / (beta * beta)); // Z^2 para proyectil (z=1 protón)
    Double_t term2 = 0.5 * TMath::Log((2 * me_c2 * beta * beta * gamma * gamma * 200.0) / (I * I)); // Simplificado
    // Nota: El '200' es un ajuste burdo de Tmax para este rango, la fórmula completa es más larga pero esto da una buena aprox.
    
    // Fórmula estándar simplificada
    Double_t dEdx_mass = (K / A) * (1/(beta*beta)) * Z * (TMath::Log(2*me*beta*beta*gamma*gamma/I) - beta*beta);
    
    // Convertir de MeV*cm2/g a MeV/cm multiplicando por densidad
    return dEdx_mass * density;
}

void CalculoBetheBloch() {
    // 1. Calcular densidad real a 300 Torr
    Double_t pressure_torr = 300.0;
    Double_t density = rho_STP * (pressure_torr / 760.0);

    std::cout << "--- PARAMETROS FISICOS ---" << std::endl;
    std::cout << "Gas: Hidrogeno (H2)" << std::endl;
    std::cout << "Presion: " << pressure_torr << " Torr" << std::endl;
    std::cout << "Densidad calculada: " << density << " g/cm3" << std::endl;
    std::cout << "--------------------------" << std::endl;

    // 2. Definir energías a probar (Protones)
    const Int_t nPoints = 5;
    Double_t energies[nPoints] = {1.0, 3.0, 5.0, 10.0, 100.0}; // MeV

    // 3. Simulación numérica (Integración por pasos)
    Double_t step_size = 0.1; // Paso de 1 mm (0.1 cm)
    
    std::cout << "\n--- RESULTADOS: RANGO DEL PROTON ---" << std::endl;
    
    for (Int_t i = 0; i < nPoints; i++) {
        Double_t E = energies[i];
        Double_t initial_E = E;
        Double_t distance = 0;

        // Bucle hasta que la energía se agote (E <= 0)
        while (E > 0.01) { // Cortamos a 10 keV
            Double_t dEdx = BetheBloch(E, density);
            
            // Si el paso es muy grande para la energía que queda, ajustamos
            if (dEdx * step_size > E) {
                distance += E / dEdx;
                E = 0;
            } else {
                E -= dEdx * step_size;
                distance += step_size;
            }
            
            // Seguridad para energías muy altas (evitar bucle infinito)
            if (distance > 20000) break; // 200 metros
        }

        std::cout << "Proton de " << initial_E << " MeV -> Recorre: " 
                  << distance << " cm (" << distance/100.0 << " metros)" << std::endl;
    }
}
