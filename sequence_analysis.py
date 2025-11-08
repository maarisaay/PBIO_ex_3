from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

def find_motif_positions(sequence, motif):
    positions = []
    start = 0
    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    return positions


def translate_all_frames(dna_seq):
    seq = Seq(dna_seq)
    reverse_seq = seq.reverse_complement()

    translations = []

    # 3 ramki odczytu dla nici forward (1, 2, 3)
    for frame in range(3):
        frame_seq = seq[frame:]
        trimmed_length = len(frame_seq) - (len(frame_seq) % 3)
        frame_seq_trimmed = frame_seq[:trimmed_length]
        protein = frame_seq_trimmed.translate(to_stop=False)
        translations.append(str(protein))

    # 3 ramki odczytu dla nici reverse (4, 5, 6)
    for frame in range(3):
        frame_seq = reverse_seq[frame:]
        trimmed_length = len(frame_seq) - (len(frame_seq) % 3)
        frame_seq_trimmed = frame_seq[:trimmed_length]
        protein = frame_seq_trimmed.translate(to_stop=False)
        translations.append(str(protein))

    return translations


def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0


def get_nucleotide_counts(sequence):
    return {
        'A': sequence.count('A'),
        'T': sequence.count('T'),
        'G': sequence.count('G'),
        'C': sequence.count('C')
    }


def create_sequence_charts(sequences, output_dir="results"):

    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "sequence_charts.png")

    # Przygotowanie danych dla pierwszych 30 pozycji
    all_sequences_30 = []
    for record in sequences:
        seq_str = str(record.seq).upper()[:30]
        all_sequences_30.append(seq_str)

    positions = range(1, 31)
    nucleotides = ['A', 'T', 'G', 'C']

    # Dane dla wykresu słupkowego - rozkład nukleotydów
    nucleotide_counts = {nuc: 0 for nuc in nucleotides}
    for seq in all_sequences_30:
        for nuc in seq:
            if nuc in nucleotide_counts:
                nucleotide_counts[nuc] += 1

    # Dane dla heatmapy - pozycja vs typ nukleotydu
    heatmap_data = np.zeros((len(nucleotides), 30))
    for seq in all_sequences_30:
        for pos, nuc in enumerate(seq):
            if nuc in nucleotides:
                nuc_idx = nucleotides.index(nuc)
                heatmap_data[nuc_idx][pos] += 1

    # Dane dla wykresu liniowego - % GC w każdej pozycji
    gc_percentages = []
    for pos in range(30):
        gc_count = 0
        total_count = 0
        for seq in all_sequences_30:
            if pos < len(seq):
                total_count += 1
                if seq[pos] in ['G', 'C']:
                    gc_count += 1
        gc_percentages.append((gc_count / total_count * 100) if total_count > 0 else 0)

    fig, axes = plt.subplots(3, 1, figsize=(12, 14))

    axes[0].bar(nucleotide_counts.keys(), nucleotide_counts.values(),
                color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
    axes[0].set_title('Rozkład nukleotydów w pierwszych 30 pozycjach', fontsize=14, fontweight='bold')
    axes[0].set_xlabel('Nukleotyd')
    axes[0].set_ylabel('Liczba wystąpień')
    axes[0].grid(axis='y', alpha=0.3)

    sns.heatmap(heatmap_data, annot=False, cmap='YlOrRd', ax=axes[1],
                xticklabels=range(1, 31), yticklabels=nucleotides, cbar_kws={'label': 'Liczba wystąpień'})
    axes[1].set_title('Heatmapa: Pozycja vs Typ nukleotydu', fontsize=14, fontweight='bold')
    axes[1].set_xlabel('Pozycja w sekwencji')
    axes[1].set_ylabel('Nukleotyd')

    axes[2].plot(range(1, 31), gc_percentages, marker='o', linewidth=2, markersize=6, color='#2ca02c')
    axes[2].set_title('Zawartość GC (%) w poszczególnych pozycjach', fontsize=14, fontweight='bold')
    axes[2].set_xlabel('Pozycja w sekwencji')
    axes[2].set_ylabel('% GC')
    axes[2].grid(True, alpha=0.3)
    axes[2].set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Wykresy zapisane do pliku: {output_file}")
    plt.close()


def create_csv_report(sequences, motifs, output_dir="results"):

    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "sequence_report.csv")

    data = []

    for record in sequences:
        seq_id = record.id
        sequence = str(record.seq).upper()
        seq_obj = Seq(sequence)

        nuc_counts = get_nucleotide_counts(sequence)

        gc_percent = calculate_gc_content(sequence)

        # Pozycje motywów
        motif_positions = {}
        for motif in motifs:
            positions = find_motif_positions(sequence, motif)
            motif_positions[f'motif_{motif}'] = str(positions)

        # Reverse complement - pierwsze 10 nukleotydów
        reverse_comp = str(seq_obj.reverse_complement())[:10]

        translations = translate_all_frames(sequence)
        protein_lengths = [len(prot) for prot in translations]

        row = {
            'sequence_id': seq_id,
            'count_A': nuc_counts['A'],
            'count_T': nuc_counts['T'],
            'count_G': nuc_counts['G'],
            'count_C': nuc_counts['C'],
            'gc_percent': round(gc_percent, 2),
            **motif_positions,
            'reverse_complement_10': reverse_comp,
            'protein_length_frame1': protein_lengths[0],
            'protein_length_frame2': protein_lengths[1],
            'protein_length_frame3': protein_lengths[2],
            'protein_length_frame4': protein_lengths[3],
            'protein_length_frame5': protein_lengths[4],
            'protein_length_frame6': protein_lengths[5]
        }
        data.append(row)

    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    print(f"\n✓ Raport CSV zapisany do pliku: {output_file}")

    return df


def calculate_balance_score(sequence):
    """
    Kryterium A: Najbardziej zbalansowana sekwencja.
    Im mniejsza różnica między nukleotydami, tym wyższy wynik.
    """
    counts = get_nucleotide_counts(sequence)
    count_values = list(counts.values())
    max_diff = max(count_values) - min(count_values)
    # Im mniejsza różnica, tym lepszy wynik (odwracamy)
    score = 1000 - max_diff
    return score


def scoring(sequence):
    """
    Kryterium B: Najciekawsza sekwencja (własna definicja).

    KRYTERIUM: "Potencjał kodujący i różnorodność genetyczna"

    Sekwencja jest oceniana jako ciekawsza jeśli:
    1. Zawiera wiele kodonów START (ATG) - potencjalne geny
    2. Ma zrównoważoną zawartość GC (40-60%) - optymalna stabilność
    3. Zawiera miejsca restrykcyjne (GAATTC) - przydatne w klonowaniu
    4. Ma niską powtarzalność (różnorodność sekwencji)
    5. Zawiera sekwencje regulatorowe (TATA box)

    Punktacja:
    - Każdy ATG: +10 punktów
    - GC w zakresie 40-60%: +50 punktów
    - Każde GAATTC: +20 punktów
    - Każda TATA: +15 punktów
    - Bonus za różnorodność: +30 punktów (jeśli wszystkie nukleotydy >15%)
    """

    score = 0

    # 1. Kodony START (ATG)
    atg_count = sequence.count('ATG')
    score += atg_count * 10

    # 2. Zrównoważona zawartość GC
    gc_content = calculate_gc_content(sequence)
    if 40 <= gc_content <= 60:
        score += 50

    # 3. Miejsca restrykcyjne (EcoRI)
    gaattc_count = sequence.count('GAATTC')
    score += gaattc_count * 20

    # 4. TATA box (sekwencje regulatorowe)
    tata_count = sequence.count('TATA')
    score += tata_count * 15

    # 5. Różnorodność nukleotydów
    counts = get_nucleotide_counts(sequence)
    total = sum(counts.values())
    percentages = [count / total * 100 for count in counts.values()]
    if all(p > 15 for p in percentages):
        score += 30

    return score


def rate_sequences(sequences, output_dir="results"):

    os.makedirs(output_dir, exist_ok=True)

    results_a = []
    results_b = []

    for record in sequences:
        seq_id = record.id
        sequence = str(record.seq).upper()

        # Kryterium A - zbalansowanie
        balance_score = calculate_balance_score(sequence)
        results_a.append((balance_score, record))

        # Kryterium B - własne
        interest_score = scoring(sequence)
        results_b.append((interest_score, record))

    results_a.sort(reverse=True, key=lambda x: x[0])
    results_b.sort(reverse=True, key=lambda x: x[0])

    top3_a = results_a[:3]
    top3_b = results_b[:3]

    print("\n" + "=" * 60)
    print("KRYTERIUM A: Najbardziej zbalansowane sekwencje (TOP 3)")
    print("=" * 60)
    for i, (score, record) in enumerate(top3_a, 1):
        counts = get_nucleotide_counts(str(record.seq).upper())
        print(f"{i}. {record.id} - Wynik: {score}")
        print(f"   Nukletydy: A={counts['A']}, T={counts['T']}, G={counts['G']}, C={counts['C']}")

    print("\n" + "=" * 60)
    print("KRYTERIUM B: Najciekawsze sekwencje (TOP 3)")
    print("Ocena oparta na potencjale kodującym i różnorodności genetycznej")
    print("=" * 60)
    for i, (score, record) in enumerate(top3_b, 1):
        seq_str = str(record.seq).upper()
        print(f"{i}. {record.id} - Wynik: {score}")
        print(f"   ATG: {seq_str.count('ATG')}, GC%: {calculate_gc_content(seq_str):.1f}, "
              f"GAATTC: {seq_str.count('GAATTC')}, TATA: {seq_str.count('TATA')}")

    file_a = os.path.join(output_dir, "kryterium_a.fasta")
    file_b = os.path.join(output_dir, "kryterium_b.fasta")

    SeqIO.write([record for _, record in top3_a], file_a, "fasta")
    SeqIO.write([record for _, record in top3_b], file_b, "fasta")

    print(f"\n✓ TOP 3 sekwencje zapisane do plików: {file_a}, {file_b}")


def analyze_sequences(fasta_file):
    """Główna funkcja analizująca sekwencje z pliku FASTA."""
    motifs = ['ATG', 'TATA', 'GAATTC']

    try:
        sequences = list(SeqIO.parse(fasta_file, "fasta"))

        if not sequences:
            print("Nie znaleziono sekwencji w pliku.")
            return

        print(f"\nWczytano {len(sequences)} sekwencji z pliku {fasta_file}")

        print("\n" + "=" * 60)
        print("PODSTAWOWA ANALIZA SEKWENCJI")
        print("=" * 60)

        for record in sequences:
            seq_id = record.id
            sequence = str(record.seq).upper()

            print(f"\nWczytano sekwencję o id {seq_id}:")
            print(sequence[:50] + "..." if len(sequence) > 50 else sequence)

            # Znajdowanie motywów
            print("\nZnalezione motywy:")
            for motif in motifs:
                positions = find_motif_positions(sequence, motif)
                print(f"{motif}:{positions}")

            # Translacja we wszystkich 6 ramkach odczytu
            print("\nTranslacja:")
            translations = translate_all_frames(sequence)
            for i, translation in enumerate(translations, 1):
                display = translation[:50] + "..." if len(translation) > 50 else translation
                print(f"Ramka {i}:{display}")

            print("-" * 60)

        print("\n" + "=" * 60)
        print("OPCJA A: TWORZENIE WYKRESÓW")
        print("=" * 60)
        create_sequence_charts(sequences)

        print("\n" + "=" * 60)
        print("OPCJA B: TWORZENIE RAPORTU CSV")
        print("=" * 60)
        create_csv_report(sequences, motifs)

        print("\n" + "=" * 60)
        print("OPCJA C: OCENA SEKWENCJI")
        print("=" * 60)
        rate_sequences(sequences)

        print("\n" + "=" * 60)
        print("ANALIZA ZAKOŃCZONA POMYŚLNIE!")
        print("=" * 60)

    except FileNotFoundError:
        print(f"Błąd: Nie znaleziono pliku '{fasta_file}'")
        print("Upewnij się, że plik znajduje się w tym samym katalogu co skrypt.")
    except Exception as e:
        print(f"Wystąpił błąd: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    fasta_file = "data/ls_orchid.fasta"

    print("=" * 60)
    print("KOMPLETNA ANALIZA SEKWENCJI DNA")
    print("=" * 60)

    analyze_sequences(fasta_file)