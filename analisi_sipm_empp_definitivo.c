//questo codice è in grado di leggere un N file di testo contenenti i dati di uno spettro e di eseguire un multifit gaussiano. 
//una volta eseguito il fit, il programma restituisce i centroidi e le sigma delle gaussiane fittate, calcola la media delle differenze
//tra i centroidi e le sigma, ed esegue un fit lineare pesato per trovare la tensione di breakdown del SiPM.
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "TH1I.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"


#define N_CHANNELS  4096

void MegaIperGaussianFit(TH1I *histogram, const std::vector<std::pair<double, double>>& params, int a, int b, std::vector<double>& centroidi, std::vector<double>& sigma) {
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

        if (isLocalMaxima && y[i]>20) {
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
        if (vettore2[i]>15 && vettore2[i]<40 && vettore2[i+1]>20 && vettore2[i+1]<40) {
            differenze_c.push_back(fabs(vettore1[i + 1] - vettore1[i]));
            differenze_s.push_back(pow(vettore2[i + 1], 2) + pow(vettore2[i], 2));
            }
    }

    double somma_c = std::accumulate(differenze_c.begin(), differenze_c.end(), 0);
    double somma_s = std::accumulate(differenze_s.begin(), differenze_s.end(), 0);
    double media_c = somma_c / differenze_c.size();
    double media_s = sqrt(somma_s) / differenze_s.size();

    return std::make_pair(media_c, media_s);
}


void analisi_sipm_empp_definitivo()
{
    for (size_t i=1; i<12; i++){
    std::ostringstream nomeFileStream;
    nomeFileStream << "Run" << i << "_PHA_HG_0_50.txt";
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
    std:: cout << "n_values: " << counts << std::endl;

    //Disegna l'istogramma
    TCanvas* c = new TCanvas();
    hspectrum->Draw(); 
    hspectrum->GetXaxis()->SetRangeUser(0, 2000); //Disegna solo l'intervallo definito
    hspectrum->GetXaxis()->SetTitle("Canali");
    hspectrum->GetYaxis()->SetTitle("Conteggi");
    hspectrum->SetLineColor(kBlue);
    hspectrum->SetFillColor(kBlue);
    hspectrum->SetFillStyle(3002);

    std::vector<std::pair<double, double>> massimiYX;

    //Trova il valore di x a partire dal quale y è maggiore o uguale a 10 (inizio intervallo)
    size_t startX = 0;
    const int threshold_m = 10;
    for (size_t i = 0; i < y.size(); ++i) {
        if (y[i] >= threshold_m) {
            startX = i;
            break;
        }
    }
    std::cout << "Inizio ricerca massimi da x = " << startX << std::endl;

    //Trova il valore di x a partire dal quale y è minore o uguale a 20 (fine intervallo)
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
    massimiYX = findLocalMaxima(x,y);

    //Stampa i risultati del massimo
    std::cout << "Massimi locali trovati:\n";
    for (const auto& maximum : massimiYX) {
        std::cout << "Coordinate (Y, X): (" << maximum.first << ", " << maximum.second << ")\n";
    }


    std::vector<double> centroidi;
    std::vector<double> sigma;

    //Funzione che esegue il multifit gaussiano che ha bisogno di un istogramma, dei massimi locali, degli estremi di ricerca e di due vettori vuoti da riempire con centroidi e sigma
    MegaIperGaussianFit(hspectrum, massimiYX, startX, endX, centroidi, sigma);

    //Stampa nuovi centroidi
    std::cout << "Contenuto del vettore centroidi: \n";
    for (const int& value : centroidi) {
        std::cout << value << "\n";
    }
    std::cout << std::endl;

    //Stampa nuovi sigma
    std::cout << "Contenuto del vettore centroidi: \n";
    for (const int& value : sigma) {
        std::cout << value << "\n";
    }
    std::cout << std::endl;

    auto risultato = calcolaMediaDifferenze(centroidi, sigma);
    std::cout << "il valore medio di distanze tra i picchi e': " << risultato.first << std::endl;
    std::cout << "l'errore relativo e': " << risultato.second << std::endl;

    c->Update();
    }

    //BEST-FIT LINEARE PESATO
    double delta_channel[5]={58.5 , 68.4444 , 79.0769 , 90.6429 , 105.615}; //(valori corrispondenti al run 7, 8, 9, 10, 11)
    double v_bias[5]={32, 32.5, 33, 33.5, 34};                              //(valori corrispondenti al run 7, 8, 9, 10, 11)
    double dc_errors[5]={21.1246, 10.3232, 9.67336, 10.8724, 13.199};       //valori ottenuti con la funzione calcolaMediaDifferenze sui sigma
    double vb_errors[5]={0.1, 0.1, 0.1, 0.1, 0.1};                          //errore di sensibilità 


    //Creo un grafico con errori
    TGraphErrors *grafico = new TGraphErrors(5, v_bias, delta_channel, vb_errors, dc_errors);

    //Definisco la funzione di best-fit (lineare: y = ax + b)
    TF1 *fitFunc = new TF1("fitFunc", "[0]*x + [1]", 30, 37);

    //Inizializzo i parametri iniziali (a e b)
    fitFunc->SetParameter(0, 15.0);
    fitFunc->SetParameter(1, -500.0);

    grafico->Fit(fitFunc, "W");

    //Estrazione del risultato del fit per i parametri "a" e "b" e relativi errori
    Double_t a = fitFunc->GetParameter(0);
    Double_t aError = fitFunc->GetParError(0);
    Double_t b = fitFunc->GetParameter(1);
    Double_t bError = fitFunc->GetParError(1);

    //Calcolo l'errore associato alla tensione di breakdown usando la teoria della propagazione degli errori
    double_t errore_vbr= (1/a)*(sqrt((pow(aError,2))*(pow(b/a, 2))+pow(bError,2)));
    cout << "Tensione di breakdown: " << -b/a << " +/- " << errore_vbr << endl;
    cout << "inverso di a: " << 1/a << endl;

    //Disegno il grafico con il fit
    TCanvas *c1 = new TCanvas("c1", "Fit Lineare Pesato", 800, 600);
    grafico->SetMarkerStyle(20);
    grafico->Draw("AP");
    

}
