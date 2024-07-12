//questo codice e' in grado di leggere N file di testo contenenti i dati di vari spettri e di eseguire un multifit gaussiano. 
//una volta eseguito il fit, il programma restituisce i centroidi, le sigma delle gaussiane fittate e le aree sootese alle singole gaussiane.
//Poi calcola il numero medio di fotoni per ogni pixel di una matrice di SiPM. Fatto cio', riempie un istogramma 2D che e' possibile
//visualizzare sia in modalita' "COLZ" sia in modalita' "LEGO2Z"
//il programma e' pensato per funzionare con i file di dati contenuti nella cartella "12_06_24_diffusori" 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "TH1I.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"

#define N_CHANNELS  4096 
#define RUN         22

void MegaIperGaussianFit(TH1I *histogram, const std::vector<std::pair<double, double>>& params, int a, int b, std::vector<double>& centroidi, std::vector<double>& sigma, std::vector<double>& areas, std::vector<double>& sigma_err) {
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
        ffit->SetParameter(i * paramsPerGaussian + 2, 29);                  //Sigma; n.b.: 29 e' un valore approssimativo che viene fuori da una prima visione dello spettro
        ffit->SetParLimits(i * paramsPerGaussian + 2, 1, 100);              //Limiti per la sigma
        ffit->SetParLimits(i * paramsPerGaussian, params[i].first, params[i].first+100); //Limiti per l'ampiezza
        ffit->SetNpx(10000);
    }
    
    //Eseguo il fit
    histogram->Fit("ffit", "R+");

    //Salvo i centroidi, le sigma, i rispettivi errori e le aree sottese alle singole gaussiane in quattro vettori distinti
    centroidi.clear();
    sigma.clear();
    sigma_err.clear();
    areas.clear();
    for (int i = 0; i < nGaussians-1; ++i) {
        centroidi.push_back(ffit->GetParameter(i * paramsPerGaussian + 1));
        sigma.push_back(ffit->GetParameter(i * paramsPerGaussian +2));
        sigma_err.push_back(ffit->GetParError(i * paramsPerGaussian +2));
        areas.push_back(fabs(sqrt(2*M_PI)*ffit->GetParameter(i * paramsPerGaussian)*ffit->GetParameter(i * paramsPerGaussian +2)));
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

void analisi_sipm_empp_h2()
{
    std::vector<double> n_entries;
    
    for (size_t i=0; i<64; i++){
    std::ostringstream nomeFileStream;
    nomeFileStream << "Run" << RUN << "_PHA_HG_0_" << i << ".txt";
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
    
    //Stampo il vettore con i conteggi per ogni pixel della matrice (per debug)
    std::cout << "n_entries: \n";
    for (const auto& entr : n_entries) {
        std::cout << entr << "\n";
    }

    std::vector<std::pair<double, double>> massimiYX;

    //Trova il valore di x a partire dal quale y e' maggiore o uguale a 10
    size_t startX = 0;
    const int threshold_m = 10;
    for (size_t i = 0; i < y.size(); ++i) {
        if (y[i] >= threshold_m) {
            startX = i;
            break;
        }
    }

    //Trova il valore di x a partire dal quale y e' minore o uguale a 20
    size_t endX = N_CHANNELS;
    const int threshold_M = 20;
    for (size_t i = x.size(); i > 0; --i) {
        if (y[i] >= threshold_M) {
            endX = i;
            break;
        }
    }

    //Cerca massimi locali in ogni intervallo e inserisce le coppie nel vettore di pair double "massimiYX"
    massimiYX = findLocalMaxima(x, y);

    //definisco i vettori che riempiro' con l'ampiezza dei centroidi, le sigma delle gaussiane e le aree sottese a ciascuna di esse
    std::vector<double> centroidi;
    std::vector<double> sigma;
    std::vector<double> areas;
    std::vector<double> sigma_err;

    //Funzione che esegue il multifit gaussiano che ha bisogno di un istogramma, dei massimi locali, degli estremi di ricerca e di quattro vettori vuoti da riempire con centroidi, sigma, errore di sigma e aree
    MegaIperGaussianFit(hspectrum, massimiYX, startX, endX, centroidi, sigma, areas, sigma_err);

    double temp=0;
    for(int l=0; l<areas.size(); l++){
        temp+=fabs(areas[l]*(l));
    }
    cout << "area per fotopicco: " << temp << endl;

    //calcolo del numero medio di fotoni per ogni pixel della matrice e riempimento del vettore per realizzaizone istogramma 2D
    n_entries.push_back(temp/counts);

    }


    //Definizione del vettore con le posizioni dei pixel fornite dal costruttore
    std::vector <int> dati_1={35, 33, 32, 34, 3, 1, 0, 2, 37, 39, 38, 36, 5, 7, 6, 4, 43, 41, 40, 42, 11, 9, 8, 10, 45, 47, 46, 44, 13, 15, 14, 12, 51, 49, 48, 50, 19, 17, 16, 18, 53, 55, 54, 52, 21, 23, 22, 20, 59, 57, 56, 58, 27, 25, 24, 26, 61, 63, 62, 60, 29, 31, 30, 28};
    
    //Creazione dell'istogramma 2D
    std::ostringstream titleStream2;
    titleStream2 << "Matrice Run" << RUN << "";
    std::string title2 = titleStream2.str();
    TH2I *h2spectrum = new TH2I("h2spectrum", title2.c_str(), 8, 0, 8, 8, 0, 8);

    //Riempimento dell'istogramma 2D tramite due cicli annidati
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            h2spectrum->Fill(i, j, n_entries[dati_1[i*8+j]]);
        }
    }

    //seguono una serie di comandi utili se si vuole cambiare la palette di colori dell'istogramma 2D
    //Definizione del range di colori (definizione palette)
    //const Int_t NRGBs = 5;
    //const Int_t NCont = 255;
    //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    //Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    //Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    //Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    //TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    //gStyle->SetNumberContours(NCont);

    gStyle->SetPalette(kRainBow);
    double min_value = *std::min_element(n_entries.begin(), n_entries.end());
    double max_value = *std::max_element(n_entries.begin(), n_entries.end());
    //std::cout << "Valore minimo in h2spectrum: " << min_value << std::endl;
    h2spectrum->GetZaxis()->SetRangeUser(0, max_value); //Imposta il range di valori di z
    h2spectrum->SetStats(false);

    //Creazione delle canvas e salvataggio su una cartella del desktop chiamata "grafici" (cambiare percorso)

    TCanvas* c1 = new TCanvas("c1", "LEGO2Z", 800, 600);
    c1->SetRightMargin(0.15);
    h2spectrum->Draw("LEGO2Z");
    std::string percorsoLEGO2Z = "/Users/mattiavalenti/Desktop/grafici/grafico_LEGO2Z_run" + std::to_string(RUN) + ".pdf";
    c1->SaveAs(percorsoLEGO2Z.c_str());
    c1->Update();

    TCanvas* c2 = new TCanvas("c2", "c2", 800, 800);
    c2->SetRightMargin(0.15);
    h2spectrum->Draw("COLZ");
    std::string percorsoCOLZ = "/Users/mattiavalenti/Desktop/grafici/grafico_COLZ_run" + std::to_string(RUN) + ".pdf";
    c2->SaveAs(percorsoCOLZ.c_str());
    c2->Update();
}
