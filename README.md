# Generator sekwencji DNA
## Funkcjonalności podstawowe:

- Generowanie losowych sekwencji DNA o dowolnej długości
- Zapis do pliku FASTA z niestandardowym ID i opisem
- Statystyki sekwencji (% nukleotydów, stosunek CG/AT)

## Dodatkowe funkcjonalności:

### Format Pro

- Podział sekwencji na linie po 60/80 znaków (standard FASTA)
- Metadata z datą generowania i nazwą narzędzia
- Eksport do dwóch formatów: FASTA lub JSON z pełnymi metadanymi

### Batch Generator

- Generowanie wielu sekwencji naraz (multi-sequence FASTA)
- Zbiorcze statystyki dla wszystkich sekwencji
- Automatyczne numerowanie sekwencji

###  Quality Control

- Walidacja zawartości GC (ostrzeżenie gdy poza zakresem 30-70%)
- Detekcja długich powtórzeń nukleotydów (>4 jednakowych z rzędu)
- Możliwość ponownej generacji przy wykryciu problemów


# Analiza sekwencji DNA

## Funkcjonalności podstawowe:

- Wczytywanie sekwencji z pliku FASTA
- Znajdowanie motywów: ATG (kodon START), TATA (TATA box), GAATTC (miejsce restrykcyjne)
-  Translacja DNA→białko w 6 ramkach odczytu (forward i reverse)

## Dodatkowe funkcjonalności:

### Sequence Charts

Generuje wykresy PNG dla pierwszych 30 pozycji sekwencji:

- Wykres słupkowy - rozkład nukleotydów
- Heatmapa - pozycja vs typ nukleotydu
- Wykres liniowy - % GC w każdej pozycji