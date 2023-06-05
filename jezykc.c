#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl-2.7.1/matrix/gsl_matrix.h>
#define MAX_WYMIAR 12

int wyznacznik_laplac(int macierz[MAX_WYMIAR][MAX_WYMIAR], int wymiar) {
    if (wymiar == 1) {
        return macierz[0][0];
    }

    int wyznacznik = 0;
    int znak = 1;
    int minor[MAX_WYMIAR][MAX_WYMIAR];

    for (int kolumna = 0; kolumna < wymiar; kolumna++) {
        int minor_wiersz = 0;
        for (int i = 1; i < wymiar; i++) {
            int minor_kolumna = 0;
            for (int j = 0; j < wymiar; j++) {
                if (j == kolumna) {
                    continue;
                }
                minor[minor_wiersz][minor_kolumna] = macierz[i][j];
                minor_kolumna++;
            }
            minor_wiersz++;
        }
        wyznacznik += znak * macierz[0][kolumna] * wyznacznik_laplac(minor, wymiar - 1);
        znak = -znak;
    }

    return wyznacznik;
}

int wyznacznik_sarrusa(int macierz[MAX_WYMIAR][MAX_WYMIAR], int wymiar) {

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
        int submacierz[MAX_WYMIAR][MAX_WYMIAR];
        int sub_wymiar = wymiar - 1;

        for (int i = 1; i < wymiar; i++) {
            int sub_kolumna = 0;
            for (int j = 0; j < wymiar; j++) {
                if (j != kolumna) {
                    submacierz[i - 1][sub_kolumna] = macierz[i][j];
                    sub_kolumna++;
                }
            }
        }

        int sub_wyznacznik = wyznacznik_sarrusa(submacierz, sub_wymiar);

        if (kolumna % 2 == 0) {
            wyznacznik += macierz[0][kolumna] * sub_wyznacznik;
        } else {
            wyznacznik -= macierz[0][kolumna] * sub_wyznacznik;
        }
    }

    return wyznacznik;
}

double wyznacznik_gsl(double macierz[MAX_WYMIAR][MAX_WYMIAR], int wymiar) {
    gsl_matrix_view mat = gsl_matrix_view_array(&macierz[0][0], wymiar, wymiar);
    gsl_permutation* perm = gsl_permutation_alloc(wymiar);
    int signum;
    double wyznacznik;

    gsl_linalg_LU_decomp(&mat.matrix, perm, &signum);
    wyznacznik = gsl_linalg_LU_det(&mat.matrix, signum);

    gsl_permutation_free(perm);

    return wyznacznik;
}

void generuj_macierze(int macierz[MAX_WYMIAR][MAX_WYMIAR], int wymiar, int wiersz, int kolumna, unsigned long long* ilosc, int min, int max, int wybor) {
    if (wiersz == wymiar) {
        (*ilosc)++;
        if (wybor == 1)
            wyznacznik_laplac(macierz, wymiar);
        if (wybor == 2)
            wyznacznik_sarrusa(macierz, wymiar);
        if (wybor == 3)
            wyznacznik_gsl(macierz, wymiar);
        return;
    }

    for (int i = min; i <= max; i++) {
        macierz[wiersz][kolumna] = i;

        if (kolumna == wymiar - 1) {
            generuj_macierze(macierz, wymiar, wiersz + 1, 0, ilosc, min, max, wybor);
        } else {
            generuj_macierze(macierz, wymiar, wiersz, kolumna + 1, ilosc, min, max, wybor);
        }
    }
}

void logika_programu(unsigned long long* ilosc) {
    int wymiar, min, max, ilosc_na_sekunde, wybor;
    int macierz[MAX_WYMIAR][MAX_WYMIAR];
    float czas;
    clock_t start, koniec;

    printf("Podaj wymiar macierzy (od 4 do 12): ");
    scanf("%d", &wymiar);
    if (wymiar < 4 || wymiar > MAX_WYMIAR) {
        printf("Nieprawidłowy wymiar macierzy.\n");
        return;
    }
    printf("Podaj zakres liczb\nmin: ");
    scanf("%d", &min);
    printf("max: ");
    scanf("%d", &max);
    printf("Jaka metoda chcesz liczyc wyznaczniki?\nMetoda Rozwiniecia Laplace'a - 1\nMetoda Sarrusa - 2\nMetoda LU(biblioteka GSL) - 3\n");
    scanf("%d", &wybor);
    if (wybor < 1 || wybor > 3) {
        printf("Nieprawidłowy wybór.\n");
        return;
    }

    printf("\nGenerowanie macierzy...\n");
    start = clock();
    generuj_macierze(macierz, wymiar, 0, 0, ilosc, min, max, wybor);
    koniec = clock();
    czas = (float)(koniec - start) / CLOCKS_PER_SEC;
    ilosc_na_sekunde = (int)(*ilosc / czas);
    printf("Czas generowania: %.3f s\n", czas);
    printf("---Wygenerowano i obliczono wyznacznik %llu macierzy co daje okolo %d macierzy na sekunde---\n", *ilosc, ilosc_na_sekunde);
}

int main() {
    int wybor;
    unsigned long long ilosc = 0;

    logika_programu(&ilosc);

    double matrix[4] = {1.0, 2.0, 3.0, 4.0};
    double determinant;

    return 0;
}
