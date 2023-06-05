#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <ctime>
#include <cfloat>
#include "Eigen/Dense"
#define MAX_DIMENSION 12


int wyznacznik_laplac(int macierz[MAX_DIMENSION][MAX_DIMENSION], int wymiar) {
    if (wymiar == 1) {
        return macierz[0][0];
    }

    int wyznacznik = 0;
    int znak = 1;
    int minor[MAX_DIMENSION][MAX_DIMENSION];

    for (int kolumna = 0; kolumna < wymiar; kolumna++) {
        int wiersz_minor = 0;
        for (int i = 1; i < wymiar; i++) {
            int kolumna_minor = 0;
            for (int j = 0; j < wymiar; j++) {
                if (j == kolumna) {
                    continue;
                }
                minor[wiersz_minor][kolumna_minor] = macierz[i][j];
                kolumna_minor++;
            }
            wiersz_minor++;
        }
        wyznacznik += znak * macierz[0][kolumna] * wyznacznik_laplac(minor, wymiar - 1);
        znak = -znak;
    }

    return wyznacznik;
}

int wyznacznik_sarrus(int macierz[MAX_DIMENSION][MAX_DIMENSION], int wymiar) {

    if (wymiar == 3) {
        int wyznacznik = 0;
        wyznacznik += macierz[0][0] * macierz[1][1] * macierz[2][2];
        wyznacznik += macierz[0][1] * macierz[1][2] * macierz[2][0];
        wyznacznik += macierz[0][2] * macierz[1][0] * macierz[2][1];
        wyznacznik -= macierz[2][0] * macierz[1][1] * macierz[0][2];
        wyznacznik -= macierz[2][1] * macierz[1][2] * macierz[0][0];
        wyznacznik -= macierz[2][2] * macierz[1][0] * macierz[0][1];
        return wyznacznik;
    }

    int wyznacznik = 0;

    for (int kolumna = 0; kolumna < wymiar; kolumna++) {
        int podmacierz[MAX_DIMENSION][MAX_DIMENSION];
        int pod_wymiar = wymiar - 1;

        for (int i = 1; i < wymiar; i++) {
            int pod_kolumna = 0;
            for (int j = 0; j < wymiar; j++) {
                if (j != kolumna) {
                    podmacierz[i - 1][pod_kolumna] = macierz[i][j];
                    pod_kolumna++;
                }
            }
        }

        int pod_wyznacznik = wyznacznik_sarrus(podmacierz, pod_wymiar);

        if (kolumna % 2 == 0) {
            wyznacznik += macierz[0][kolumna] * pod_wyznacznik;
        } else {
            wyznacznik -= macierz[0][kolumna] * pod_wyznacznik;
        }
    }

    return wyznacznik;
}

int wyznacznik_eigen(int macierz[MAX_DIMENSION][MAX_DIMENSION], int wiersze, int kolumny) {
  Eigen::MatrixXd macierz_eigen(wiersze, kolumny);
  
  for (int i = 0; i < wiersze; ++i) {
    for (int j = 0; j < kolumny; ++j) {
      macierz_eigen(i, j) = static_cast<double>(macierz[i][j]);
    }
  }
  
  int wyznacznik = static_cast<int>(macierz_eigen.determinant());
  return wyznacznik;
}

void wypisz_macierz(int macierz[MAX_DIMENSION][MAX_DIMENSION], int wymiar, std::fstream& plik) {
    for (int i = 0; i < wymiar; i++) {
        for (int j = 0; j < wymiar; j++) {
            plik << macierz[i][j] << " ";
        }
        plik << "\n";
    }
    plik << "\n";
}

void oblicz_wyznaczniki_lap(int macierz[MAX_DIMENSION][MAX_DIMENSION], int wymiar, std::vector<int>& wyznaczniki) {
    wyznaczniki.push_back(wyznacznik_laplac(macierz, wymiar));
}

void oblicz_wyznaczniki_sarr(int macierz[MAX_DIMENSION][MAX_DIMENSION], int wymiar, std::vector<int>& wyznaczniki) {
    wyznaczniki.push_back(wyznacznik_sarrus(macierz, wymiar));
}

void oblicz_wyznaczniki_eigen(int macierz[MAX_DIMENSION][MAX_DIMENSION], int wymiar, std::vector<int>& wyznaczniki) {
    wyznaczniki.push_back(wyznacznik_eigen(macierz, wymiar, wymiar));
}

void generuj_i_wypisz_macierze(int macierz[MAX_DIMENSION][MAX_DIMENSION], int wymiar, int wiersz, int kolumna, unsigned long long& ilosc, std::fstream& plik, std::vector<int>& wyznaczniki, int min, int max) {
    if (wiersz == wymiar) {
        ilosc += 1;
        oblicz_wyznaczniki_sarr(macierz, wymiar, wyznaczniki);
        wypisz_macierz(macierz, wymiar, plik);
        return;
    }

    for (int i = min; i <= max; i++) {
        macierz[wiersz][kolumna] = i;

        if (kolumna == wymiar - 1) {
            generuj_i_wypisz_macierze(macierz, wymiar, wiersz + 1, 0, ilosc, plik, wyznaczniki, min, max);
        } else {
            generuj_i_wypisz_macierze(macierz, wymiar, wiersz, kolumna + 1, ilosc, plik, wyznaczniki, min, max);
        }
    }
}

void generuj_macierze(int macierz[MAX_DIMENSION][MAX_DIMENSION], int wymiar, int wiersz, int kolumna, unsigned long long& ilosc, int min, int max, int wybor, std::vector<int>& wyznaczniki) {
    if (wiersz == wymiar) {
        ilosc += 1;
        if (wybor == 1)
            oblicz_wyznaczniki_lap(macierz, wymiar, wyznaczniki);
        if (wybor == 2)
            oblicz_wyznaczniki_sarr(macierz, wymiar, wyznaczniki);
        if(wybor == 3)
            oblicz_wyznaczniki_eigen(macierz, wymiar, wyznaczniki);
        return;
    }

    for (int i = min; i <= max; i++) {
        macierz[wiersz][kolumna] = i;

        if (kolumna == wymiar - 1) {
            generuj_macierze(macierz, wymiar, wiersz + 1, 0, ilosc, min, max, wybor, wyznaczniki);
        } else {
            generuj_macierze(macierz, wymiar, wiersz, kolumna + 1, ilosc, min, max, wybor, wyznaczniki);
        }
    }
}

void generuj_histogram(std::vector<int>& wyznaczniki, const std::string& plik) {
    std::cout << "\nGenerowanie histogramu...\n";
    std::map<int, int> histogram;
    for (int num : wyznaczniki) {
        histogram[num]++;
    }

    std::ofstream file(plik);

    for (const auto& entry : histogram) {
        file << entry.first << " " << entry.second << std::endl;
    }

    file.close();

    std::cout << "Histogram zostal zapisany do pliku: " << plik <<"\n\n";
}

void generuj_wykres(const std::string& plik, int wymiar){
    
    std::ofstream komendy("histogram-graficznie.gp");
    komendy << "set style data histogram\n";
    komendy << "set style fill solid\n";
    komendy << "set title \"Histogram wyznacznikow macierzy " << wymiar << " x " << wymiar << "\"\n";
    komendy << "plot \"" << plik << "\" using 2:xtic(1) with histogram\n";
    komendy << "pause mouse close\n";
    komendy.close();
}

void logika_programu(unsigned long long& ilosc, std::vector <int>& wyznaczniki) {
    int wymiar, min, max, ilosc_na_sekunde, wybor;
    std::fstream file;
    int macierz[MAX_DIMENSION][MAX_DIMENSION] = {};
    float czas;
    std::clock_t start, koniec;

    std::cout << "Podaj wymiar macierzy(od 4 do 12): ";
    std::cin >> wymiar;
    if (wymiar < 4 || wymiar > MAX_DIMENSION) {
        std::cout << "Nieprawidlowy wymiar macierzy.\n";
        return;
    }
    std::cout << "Podaj zakres liczb\nmin: ";
    std::cin >> min;
    std::cout << "max: ";
    std::cin >> max;
    std::cout << "Jaka metoda chcesz liczyc wyznaczniki?\nMetoda Rozwiniecia Laplacea - 1\nMetoda Sarrusa - 2\nMetoda Gaussa-Jordana(biblioteka eigen) - 3\n";
    std::cin >> wybor;
    if (wybor < 1 || wybor > 3) {
        std::cout << "Nieprawidlowy wybor.\n";
        return;
    }

    std::cout << "\nGenerowanie macierzy...\n";
    start = std::clock();
    
    generuj_macierze(macierz, wymiar, 0, 0, ilosc, min, max, wybor, wyznaczniki);
    koniec = std::clock();
    czas = (koniec - start) / (double)CLOCKS_PER_SEC;
    ilosc_na_sekunde = ilosc/czas;
    std::cout << "Czas generowania: " << czas << " s\n";
    std::cout<<"---Wygenerowano i obliczono wyznacznik "<<ilosc<<" macierzy co daje okolo " << ilosc_na_sekunde <<" macierzy na sekunde---\n";
 
    std::cout << "Czy chcesz wypisac macierze do pliku i stworzyc histogram? (rekomendowane tylko jesli wygenerowano mniej niz 100 000 000 macierzy)\nTAK - 1    NIE - 2\n";
    std::cin >> wybor;
    if (wybor < 1 || wybor > 2) {
        std::cout << "Nieprawidlowy wybor.\n";
        return;
    }

    if(wybor == 1){
        file.open("macierze.txt", std::ios::out);
        std::cout << "\nWypisywanie macierzy...\n";
        generuj_i_wypisz_macierze(macierz, wymiar, 0, 0, ilosc, file, wyznaczniki, min, max);
        std::cout << "Macierze zostaly zapisane do pliku: macierze.txt\n";
        generuj_histogram(wyznaczniki, "histogram.txt");
        generuj_wykres("histogram.txt", wymiar);
        file.close();
    }
}

int main() {
    int wybor;
    unsigned long long ilosc = 0;
    std::vector <int> wyznaczniki;

    logika_programu(ilosc, wyznaczniki);

    return 0;
}
