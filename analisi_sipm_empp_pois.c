//questo codice è in grado di leggere un file di testo contenente i dati di uno spettro e di eseguire un multifit gaussiano. 
//una volta eseguito il fit, il programma restituisce i centroidi e le sigma delle gaussiane fittate, calcola la media delle differenze
//tra i centroidi e le sigma, ed esegue un fit poissoniano per verificare che i dati sperimentali seguano la distribuzione di Poisson.
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "TH1I.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"

#define N_CHANNELS  4096

void MegaIperGaussianFit(TH1I *histogram, const std::vector<std::pair<double, double>>& params, int a, int b, std::vector<double>& centroidi, std::vector<double>& sigma, std::vector<double>& areas) {
    const int paramsPerGaussian = 3;  //Parametri per ogni gaussiana (ampiezza, centroide, sigma) 
    int nGaussians=params.size();    
    
    //Creo la stringa della funzione di fit dinamicamente
    std::string formula;
    for (int i = 0; i < nGaussians; ++i) {
        if (i > 0) {
            formula += "+";
        }
        formula += "gaus(" + std::to_string(i * paramsPerGaussian) + ")";
    }

    //Creo la funzione di fit
    TF1 *ffit = new TF1("ffit", formula.c_str(), a, b); //gli ultimi numeri servono per definire gli estremi dell'intervallo in cui effettuare il fit

    //Imposto i parametri iniziali
    for (int i = 0; i < nGaussians; ++i) {
        ffit->SetParameter(i * paramsPerGaussian, params[i].first);         //Ampiezza
        ffit->SetParameter(i * paramsPerGaussian + 1, params[i].second);    //Centroide
        ffit->SetParameter(i * paramsPerGaussian + 2, 29);                  //Sigma; n.b.: 29 è un valore approssimativo che viene fuori da una prima visione dello spettro
        ffit->SetNpx(10000);
    }
    
    //Eseguo il fit
    histogram->Fit("ffit", "R+");

    //Salvo i centroidi e le sigma in due vettori distinti
    centroidi.clear();
    sigma.clear();
    for (int i = 0; i < nGaussians-1; ++i) {
        centroidi.push_back(ffit->GetParameter(i * paramsPerGaussian + 1));
        sigma.push_back(ffit->GetParameter(i * paramsPerGaussian +2));
        areas.push_back(sqrt(2*M_PI)*ffit->GetParameter(i * paramsPerGaussian)*ffit->GetParameter(i * paramsPerGaussian +2));
    }
}

std::vector<std::pair<double, double>> findLocalMaxima(const std::vector<double>& x, const std::vector<double>& y, size_t range = 23) {
    if (x.size() != y.size() || x.size() == 0) {
        throw std::invalid_argument("Invalid input sizes.");
    }
    
    std::vector<std::pair<double, double>> localMaxima;
    int n_max = 0;

    for (size_t i = 0; i < y.size(); ++i) {
        bool isLocalMaxima = true;
        
        size_t start = (i < range) ? 0 : i - range;
        size_t end = (i + range >= y.size()) ? y.size() - 1 : i + range;

        for (size_t j = start; j <= end; ++j) {
            if (j != i && y[j] >= y[i]) {
                isLocalMaxima = false;
                break;
            }
        }

        if (isLocalMaxima && y[i] > 20) {
            localMaxima.push_back(std::make_pair(y[i], x[i]));
            n_max++;
        }    
    }

    return localMaxima;
}

std::pair<double, double> calcolaMediaDifferenze(const std::vector<double>& vettore1, const std::vector<double>& vettore2) {
    if (vettore1.size() != vettore2.size() || vettore1.size() < 2) {
        std::cerr << "I vettori devono contenere almeno due elementi." << std::endl;
        return std::make_pair(0, 0);
    }

    std::vector<double> differenze_c;
    std::vector<double> differenze_s;
    for (size_t i = 0; i < vettore1.size() - 1; ++i) {
        if (vettore2[i] > 15 && vettore2[i] < 40 && vettore2[i + 1] > 20 && vettore2[i + 1] < 40) {
            differenze_c.push_back(fabs(vettore1[i + 1] - vettore1[i]));
            differenze_s.push_back(pow(vettore2[i + 1], 2) + pow(vettore2[i], 2));
        }
    }

    double somma_c = std::accumulate(differenze_c.begin(), differenze_c.end(), 0.0);
    double somma_s = std::accumulate(differenze_s.begin(), differenze_s.end(), 0.0);
    double media_c = somma_c / differenze_c.size();
    double media_s = sqrt(somma_s) / differenze_s.size();

    return std::make_pair(media_c, media_s);
}

void analisi_sipm_empp_pois() {
    std::string nomeFile = "Run8_PHA_HG_0_50.txt";
    
    //Apertura del file di input
    std::ifstream inputFile(nomeFile);
    if (!inputFile.is_open()) {
        std::cerr << "Errore nell'apertura del file: " << nomeFile << std::endl;
        return;
    }
    
    //Creazione dell'istogramma
    TH1I *hspectrum = new TH1I("hspectrum", "Spettro Run 8 - Illuminazione 6", 4096, 0, 4095);

    int value;
    int counts = 0;
    std::vector<double> x;
    std::vector<double> y;

    for (int i = 0; i < N_CHANNELS; i++) {
        if (!(inputFile >> value)) {
            std::cerr << "Errore nella lettura del file alla riga " << i + 1 << std::endl;
        }
        hspectrum->SetBinContent(i + 1, value);
        counts += value;
        x.push_back(i);
        y.push_back(value);
    }
    std::cout << "n_values: " << counts << std::endl;

    //Disegna l'istogramma
    TCanvas* c = new TCanvas();
    hspectrum->Draw();
    hspectrum->GetXaxis()->SetRangeUser(0, 2000);
    hspectrum->GetXaxis()->SetTitle("Canali");
    hspectrum->GetYaxis()->SetTitle("Conteggi");
    hspectrum->SetLineColor(kBlue);
    hspectrum->SetFillColor(kBlue);
    hspectrum->SetFillStyle(3002);

    std::vector<std::pair<double, double>> massimiYX;

    //Trova il valore di x a partire dal quale y è maggiore o uguale a 10
    size_t startX = 0;
    const int threshold_m = 10;
    for (size_t i = 0; i < y.size(); ++i) {
        if (y[i] >= threshold_m) {
            startX = i;
            break;
        }
    }
    std::cout << "Inizio ricerca massimi da x = " << startX << std::endl;

    //Trova il valore di x a partire dal quale y è minore o uguale a 20
    size_t endX = N_CHANNELS;
    const int threshold_M = 20;
    for (size_t i = x.size(); i > 0; --i) {
        if (y[i] >= threshold_M) {
            endX = i;
            break;
        }
    }
    std::cout << "Finisco la ricerca dei massimi da x = " << endX << std::endl;

    //Cerca massimi locali in ogni intervallo e inserisce le coppie nel vettore di pair double "massimiYX"
    massimiYX = findLocalMaxima(x, y);

    //Stampa i risultati del massimo
    std::cout << "Massimi locali trovati:\n";
    for (const auto& maximum : massimiYX) {
        std::cout << "Coordinate (Y, X): (" << maximum.first << ", " << maximum.second << ")\n";
    }
    
    //definisco i vettori che riempirò con l'ampiezza dei centroidi, le sigma delle gaussiane e le aree sottese a ciascuna di esse
    std::vector<double> centroidi;
    std::vector<double> sigma;
    std::vector<double> areas;

    //Funzione che esegue il multifit gaussiano che ha bisogno di un istogramma, dei massimi locali, degli estremi di ricerca e di due vettori vuoti da riempire con centroidi e sigma
    MegaIperGaussianFit(hspectrum, massimiYX, startX, endX, centroidi, sigma, areas);

    //Stampa nuovi centroidi
    std::cout << "Contenuto del vettore centroidi: \n";
    for (const auto& value : centroidi) {
        std::cout << value << "\n";
    }

    //Stampa nuovi sigma
    std::cout << "Contenuto del vettore sigma: \n";
    for (const auto& value : sigma) {
        std::cout << value << "\n";
    }

    //Stampa le aree
    std::cout << "Aree delle gaussiane: \n";
    for (const auto& area : areas) {
        std::cout << area << "\n";
    }

    auto risultato = calcolaMediaDifferenze(centroidi, sigma);
    std::cout << "Il valore medio di distanze tra i picchi è: " << risultato.first << std::endl;
    std::cout << "L'errore relativo è: " << risultato.second << std::endl;

    //Definisco il vettore dei fotoelettroni e il vettore degli errori delle aree
    double photo_el[16]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
    std::vector<double> error_areas;

    //Riempio il vettore degli errori delle aree assumendo che l'errore sia la radice di 2*pi*sigma (si assume che delta_ampiezza=1)
    for(int i=0; i<areas.size(); i++){
        error_areas.push_back(sqrt(2*M_PI)*sigma[i]);
    }

    //Stampa l'errore delle aree
    std::cout << "errori delle aree: \n";
    for (const auto& area : error_areas) {
        std::cout << area << "\n";
    }
    
    //Fitto i punti con una funzione di Poisson
    //il "-1" è dovuto all'esclusione dei valori cossrispondeti all'ultimo picco, il cui fit non quadra mai (spero di risolvere)
    TGraphErrors *graph = new TGraphErrors(areas.size()-1, photo_el, areas.data()-1, 0, error_areas.data()-1);
    
    //Creo la funzione di fit
    TF1 *f1 = new TF1("f1","[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", 0, 15); // "xmin" = 0, "xmax" = 15
    f1->SetParameters(1, 1, 1); //Non si possono settare valori iniziali nulli!!
    f1->SetLineColor(kAzure+4);
    graph->Fit("f1", "R"); // "R" = fit tra "xmin" e "xmax" della funzione "f1"

    TCanvas *c2 = new TCanvas("c2", "Fit aree gaussiane", 800, 600);
    c2->cd();
    graph->SetMarkerStyle(83);
    graph->SetMarkerColor(kViolet+4);
    graph->GetXaxis()->SetTitle("Numero di fotoelettroni");
    graph->GetYaxis()->SetTitle("Conteggi (area)");
    graph->Draw("AP");
}
