
import csv
import matplotlib.pyplot as plt
from pathlib import Path
import tkinter as tk
from tkinter import filedialog


# ==========================================
# GC CONTENT ANALYZER
# ==========================================

def get_fasta_file():

    root = tk.Tk()
    root.withdraw()

    filename = filedialog.askopenfilename(
        title="Select FASTA File",
        filetypes=[
            ("FASTA files", "*.fasta *.fa *.fna *.fas"),
            ("All files", "*.*")
        ]
    )

    root.destroy()

    if not filename:
        return None

    with open(filename, "r") as file:
        return file.read()


# ==========================================
# PASTE FASTA SEQUENCE
# ==========================================

def paste_fasta():

    print("\nPaste your FASTA sequence below.")
    print("Type END on a separate line when finished.\n")

    lines = []

    while True:

        line = input()

        if line.strip().upper() == "END":
            break

        lines.append(line)

    return "\n".join(lines)


# ==========================================
# PARSE FASTA
# ==========================================

def parse_fasta(text):

    sequences = {}
    current_name = None

    for line in text.splitlines():

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            current_name = line[1:].strip()

            if not current_name:
                raise ValueError("Empty FASTA header.")

            if current_name in sequences:
                raise ValueError(
                    f"Duplicate sequence name: {current_name}"
                )

            sequences[current_name] = ""

        elif current_name is None:

            raise ValueError(
                "FASTA header missing. Each sequence must start with >"
            )

        else:

            sequences[current_name] += line.upper()

    if not sequences:
        raise ValueError("No FASTA sequences found.")

    return sequences


# ==========================================
# DNA SEQUENCE ANALYSIS
# ==========================================

def analyze_sequence(name, sequence):

    sequence = "".join(sequence.split()).upper()

    counts = {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C")
    }

    length = len(sequence)

    valid_bases = sum(counts.values())

    gc_count = counts["G"] + counts["C"]

    at_count = counts["A"] + counts["T"]

    if valid_bases > 0:

        gc_content = (gc_count / valid_bases) * 100

        at_content = (at_count / valid_bases) * 100

    else:

        gc_content = 0
        at_content = 0

    invalid_count = length - valid_bases

    return {
        "Sequence": name,
        "Length": length,
        "A": counts["A"],
        "T": counts["T"],
        "G": counts["G"],
        "C": counts["C"],
        "GC_Content_%": round(gc_content, 2),
        "AT_Content_%": round(at_content, 2),
        "Invalid_or_Ambiguous": invalid_count
    }


# ==========================================
# SAVE CSV REPORT
# ==========================================

def save_report(results):

    filename = "GC_Analysis_Report.csv"

    with open(filename, "w", newline="") as file:

        writer = csv.DictWriter(
            file,
            fieldnames=results[0].keys()
        )

        writer.writeheader()

        writer.writerows(results)

    print("\nCSV report saved:", filename)


# ==========================================
# GENERATE GRAPH
# ==========================================

def generate_graph(results):

    names = [r["Sequence"] for r in results]

    gc_values = [r["GC_Content_%"] for r in results]

    plt.figure(figsize=(10, 5))

    plt.bar(names, gc_values)

    plt.xlabel("DNA Sequences")

    plt.ylabel("GC Content (%)")

    plt.title("GC Content Analysis")

    plt.ylim(0, 100)

    plt.xticks(rotation=45, ha="right")

    plt.tight_layout()

    plt.savefig("GC_Content_Graph.png", dpi=300)

    plt.show()

    print("Graph saved: GC_Content_Graph.png")


# ==========================================
# DISPLAY RESULTS
# ==========================================

def display_results(results):

    print("\n" + "=" * 55)

    print("           DNA ANALYSIS RESULTS")

    print("=" * 55)

    for result in results:

        print("\nSequence:", result["Sequence"])

        print("Length:", result["Length"], "bp")

        print("A:", result["A"])

        print("T:", result["T"])

        print("G:", result["G"])

        print("C:", result["C"])

        print("GC Content:",
              result["GC_Content_%"], "%")

        print("AT Content:",
              result["AT_Content_%"], "%")

        print("Invalid/Ambiguous:",
              result["Invalid_or_Ambiguous"])

        print("-" * 40)


# ==========================================
# MAIN PROGRAM
# ==========================================

def main():

    print("\n" + "=" * 50)

    print("          GC CONTENT ANALYZER")

    print("=" * 50)

    print("\nChoose your input method:")

    print("1. Upload FASTA file")

    print("2. Copy and paste FASTA sequence")

    choice = input("\nEnter your choice (1 or 2): ").strip()

    if choice == "1":

        fasta_text = get_fasta_file()

        if fasta_text is None:

            print("No file selected. Exiting.")

            return

    elif choice == "2":

        fasta_text = paste_fasta()

    else:

        print("Invalid choice. Please select 1 or 2.")

        return

    try:

        sequences = parse_fasta(fasta_text)

        results = []

        for name, sequence in sequences.items():

            results.append(
                analyze_sequence(name, sequence)
            )

        display_results(results)

        save_report(results)

        generate_graph(results)

        print("\nAnalysis completed successfully!")

    except (ValueError, OSError) as error:

        print("\nError:", error)


# ==========================================
# RUN PROGRAM
# ==========================================

if __name__ == "__main__":

    main()
