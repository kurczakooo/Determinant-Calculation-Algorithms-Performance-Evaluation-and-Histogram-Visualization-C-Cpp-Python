set style data histogram
set style fill solid
set title "Histogram wyznacznikow macierzy 4 x 4"
plot "histogram.txt" using 2:xtic(1) with histogram
pause mouse close
