import time
import numpy

MAX_WYMIAR = 12

def wyznacznik_laplac(macierz, wymiar):
    if wymiar == 1:
        return macierz[0][0]
    
    wyznacznik = 0
    znak = 1
    macierz_dopelnien = [[0] * MAX_WYMIAR for _ in range(MAX_WYMIAR)]

    for kolumna in range(wymiar):
        wiersz_dopelnien = 0
        for i in range(1, wymiar):
            kolumna_dopelnien = 0
            for j in range(wymiar):
                if j == kolumna:
                    continue
                macierz_dopelnien[wiersz_dopelnien][kolumna_dopelnien] = macierz[i][j]
                kolumna_dopelnien += 1
            wiersz_dopelnien += 1
        
        wyznacznik += znak * macierz[0][kolumna] * wyznacznik_laplac(macierz_dopelnien, wymiar - 1)
        znak = -znak
    
    return wyznacznik

def wyznacznik_sarrus(macierz, wymiar):

    if wymiar == 3:
        wyznacznik = 0
        wyznacznik += macierz[0][0] * macierz[1][1] * macierz[2][2]
        wyznacznik += macierz[0][1] * macierz[1][2] * macierz[2][0]
        wyznacznik += macierz[0][2] * macierz[1][0] * macierz[2][1]
        wyznacznik -= macierz[2][0] * macierz[1][1] * macierz[0][2]
        wyznacznik -= macierz[2][1] * macierz[1][2] * macierz[0][0]
        wyznacznik -= macierz[2][2] * macierz[1][0] * macierz[0][1]
        return wyznacznik
    
    wyznacznik = 0

    for kolumna in range(wymiar):
        podmacierz = [[0] * MAX_WYMIAR for _ in range(MAX_WYMIAR)]
        podwymiar = wymiar - 1

        for i in range(1, wymiar):
            kolumna_podmacierzy = 0
            for j in range(wymiar):
                if j != kolumna:
                    podmacierz[i - 1][kolumna_podmacierzy] = macierz[i][j]
                    kolumna_podmacierzy += 1
        
        wyznacznik_podmacierzy = wyznacznik_sarrus(podmacierz, podwymiar)

        if kolumna % 2 == 0:
            wyznacznik += macierz[0][kolumna] * wyznacznik_podmacierzy
        else:
            wyznacznik -= macierz[0][kolumna] * wyznacznik_podmacierzy
    
    return wyznacznik

def wyznacznik_numpy(macierz):  
    wyznacznik = numpy.linalg.det(macierz);
    return wyznacznik

def generuj_macierze(macierz, wymiar, wiersz, kolumna, ilosc, min_val, max_val, wybor):
    if wiersz == wymiar:
        ilosc[0] += 1
        if wybor == 1:
            wyznacznik_laplac(macierz, wymiar)
        if wybor == 2:
            wyznacznik_sarrus(macierz, wymiar)
        if wybor == 3:
            wyznacznik_numpy(macierz)
        return
    
    for i in range(min_val, max_val + 1):
        macierz[wiersz][kolumna] = i

        if kolumna == wymiar - 1:
            generuj_macierze(macierz, wymiar, wiersz + 1, 0, ilosc, min_val, max_val, wybor)
        else:
            generuj_macierze(macierz, wymiar, wiersz, kolumna + 1, ilosc, min_val, max_val, wybor)

def logika_programu():
    wymiar = int(input("Podaj wymiar macierzy (od 4 do 12): "))
    if wymiar < 4 or wymiar > MAX_WYMIAR:
        print("Nieprawidłowy wymiar macierzy.")
        return
    min_val = int(input("Podaj zakres liczb\nmin: "))
    max_val = int(input("max: "))
    wybor = int(input("Jaka metoda chcesz liczyc wyznaczniki?\nMetoda Rozwinięcia Laplace'a - 1\nMetoda Sarrusa - 2\nMetoda Gaussa Jordana(numpy) - 3\n"))
    if wybor < 1 or wybor > 3:
        print("Nieprawidłowy wybór.")
        return

    macierz = [[0] * MAX_WYMIAR for _ in range(MAX_WYMIAR)]
    ilosc = [0]
    
    print("\nGenerowanie macierzy...")
    start = time.time()
    generuj_macierze(macierz, wymiar, 0, 0, ilosc, min_val, max_val, wybor)
    koniec = time.time()
    czas = koniec - start
    ilosc_na_sekunde = int(ilosc[0] / czas)
    print("Czas generowania: %.3f s" % czas)
    print("---Wygenerowano i obliczono wyznacznik %d macierzy, co daje okolo %d macierzy na sekunde---" % (ilosc[0], ilosc_na_sekunde))

logika_programu()
