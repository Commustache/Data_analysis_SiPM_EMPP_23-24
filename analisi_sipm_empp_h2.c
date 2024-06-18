#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "TH1I.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"

#define N_CHANNELS  4096

void analisi_sipm_empp_h2()
{
    std::vector<double> n_entries;
    
    for (size_t i=0; i<64; i++){
    std::ostringstream nomeFileStream;
    nomeFileStream << "Run2_PHA_HG_0_" << i << ".txt";
    std::string nomeFile = nomeFileStream.str();
    
    //Apertura del file di input
    std::ifstream inputFile(nomeFile);
    if (!inputFile.is_open()) {
        std::cerr << "Errore nell'apertura del file: " << nomeFile << std::endl;
        continue;
    }
    
    //Creazione dell'istogramma
    std::ostringstream titleStream;
    titleStream << "Spettro Run" << i << " - Illuminazione = 6";
    std::string title = titleStream.str();
    TH1I *hspectrum = new TH1I("hspectrum", title.c_str(), 4096, 0, 4095);

    int value;
    int counts=0;
    //Inizializzazione vettori per trovare il minimo (e/o massimo) locale
    std::vector<double> x;
    std::vector<double> y;

    for (int i = 0; i < N_CHANNELS; i++) {
        if (!(inputFile >> value)) {
            std::cerr << "Errore nella lettura del file alla riga " << i+1 << std::endl;
        }
        hspectrum->SetBinContent(i+1, value); //Riempie il bin corrispondente al canale con il valore letto
        counts = counts + value;
        x.push_back(i);
        y.push_back(value);
    }
    std:: cout << "n_values_" << i << ": " << counts << std::endl;
    n_entries.push_back(counts);
    

    //Disegna l'istogramma
    //TCanvas* c = new TCanvas();
    //hspectrum->Draw(); 
    //hspectrum->GetXaxis()->SetRangeUser(0, 800); //Disegna solo l'intervallo definito
    //hspectrum->GetXaxis()->SetTitle("Canali");
    //hspectrum->GetYaxis()->SetTitle("Conteggi");
    //hspectrum->SetLineColor(kBlue);
    //hspectrum->SetFillColor(kBlue);
    //hspectrum->SetFillStyle(3002);
    
    //c->Update();
    }

    //Stampo il vettore con i conteggi per ogni pixel della matrice
    std::cout << "n_entries: \n";
    for (const auto& area : n_entries) {
        std::cout << area << "\n";
    }

    TCanvas* c = new TCanvas();

    //Definizione del vettore con le posizioni dei pixel fornite dal costruttore
    std::vector <int> dati_1={35, 33, 32, 34, 3, 1, 0, 2, 37, 39, 38, 36, 5, 7, 6, 4, 43, 41, 40, 42, 11, 9, 8, 10, 45, 47, 46, 44, 13, 15, 14, 12, 51, 49, 48, 50, 19, 17, 16, 18, 53, 55, 54, 52, 21, 23, 22, 20, 59, 57, 56, 58, 27, 25, 24, 26, 61, 63, 62, 60, 29, 31, 30, 28};
    
    //Creazione dell'istogramma 2D
    TH2I *h2spectrum = new TH2I("h2spectrum", "Spettro Run 1 - Illuminazione 6", 8, 0, 8, 8, 0, 8);

    //Riempimento dell'istogramma 2D tramite due cicli annidati
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            h2spectrum->Fill(i, j, n_entries[dati_1[i*8+j]]);
        }
    }

    //Definizione del range di colori (definizione palette)
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    h2spectrum->GetZaxis()->SetRangeUser(998800, 1000000); //Imposta il range di valori di z
    h2spectrum->Draw("COLZ");
}
