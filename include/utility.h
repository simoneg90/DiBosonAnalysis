//Questo programma contiene funzioni utili all'esecuzione di analisi

#include <stdlib.h>
#include <iostream>
#include <string.h>
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TCanvas.h"

//Questa funzione serve per creare spazio tra una stampa ed un'altra
void breakLine();
//Questa funzione stampa una stringa all'interno di un riquadro
void frame(std::string stringa);
//Questa funzione implementa le caratteristiche di un istogramma TH1D
void TH1D_config(TH1D *histo, const char* title, const char* title_x, const char* title_y, int color);
//Questa funzione implementa le caratteristiche di un istogramma TH1F
void TH1F_config(TH1F *histo, const char* title, const char* title_x, const char* title_y, int color);
//Questa funzione implementa le caratteristiche di un istogramma TH2D
void TH2D_config(TH2D *histo, const char* title, const char* title_x, const char* title_y, int color);
//Questa funzione implementa le caratteristiche di un istogramma TH2F
void TH2F_config(TH2F *histo, const char* title, const char* title_x, const char* title_y, int color);
//Questa funzione analizza gli istogrammi all'interno di un file.root e li scrive su un file di output
void LS_config(std::string rootFile, std::string outfile);
