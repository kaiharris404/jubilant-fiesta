import os
import gzip
from Bio import SeqIO
import matplotlib.pyplot as plt

# ==========================================
# MODULE 1: Parse Data (Handles .gz and .fastq)
# ==========================================
def analyze_fastq(filename):
    """Reads FASTQ data and returns sequence counts, lengths, and N-content."""
    seq_counts = {}
    lengths = []
    n_percents = []

    if filename.endswith(".gz"):
        handle = gzip.open(filename, "rt")
    else:
        handle = open(filename, "r")

    for record in SeqIO.parse(handle, "fastq"):
        seq = str(record.seq).upper()
        
        # Abundance counting (Dictionary)
        seq_counts[seq] = seq_counts.get(seq, 0) + 1
            
        # Metrics for lists
        lengths.append(len(seq))
        n_percent = (seq.count('N') / len(seq)) * 100
        n_percents.append(n_percent)

    handle.close() 
    return seq_counts, lengths, n_percents

# ==========================================
# MODULE 2: Sorting & Saving
# ==========================================
def save_top_sequences(sorted_seqs, output_filename, total_reads, top_n=15):
    """Saves the analysis results to a text file."""
    with open(output_filename, 'w') as file:
        file.write("BIOINFORMATICS ANALYSIS REPORT\n")
        file.write("=" * 30 + "\n")
        file.write(f"Total Reads Processed: {total_reads}\n")
        file.write(f"Unique Sequences Found: {len(sorted_seqs)}\n\n")
        file.write(f"Top {top_n} Sequences:\n")
        for seq, count in sorted_seqs[:top_n]:
            file.write(f"Count: {count} | Seq: {seq}\n")

def create_plots(lengths, n_percents, folder_name):
    """Generates and shows histograms."""
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.hist(lengths, bins=30, color='skyblue', edgecolor='black')
    plt.title(f'Length Distribution: {folder_name}')
    plt.xlabel('Length (bp)')
    plt.ylabel('Frequency')
    
    plt.subplot(1, 2, 2)
    plt.hist(n_percents, bins=30, color='salmon', edgecolor='black')
    plt.title(f'N-Content: {folder_name}')
    plt.xlabel('% of N bases')
    plt.ylabel('Frequency')
    
    plt.tight_layout()
    plt.show()

# ==========================================
# MAIN EXECUTION BLOCK (Interactive Mode)
# ==========================================
print("--- Bioinformatics Folder Analyzer ---")

# This asks YOU for the folder name in the terminal
target_folder = input("Enter the folder name to analyze (e.g., barcode09): ").strip()

# Check if the folder actually exists
if os.path.exists(target_folder) and os.path.isdir(target_folder):
    report_name = f"{target_folder}_analysis_report.txt"
    
    final_counts = {}
    final_lengths = []
    final_n_percents = []
    total_reads_count = 0

    print(f"\nProcessing files inside '{target_folder}'...")

    # Loop through the files
    for file_name in os.listdir(target_folder):
        if file_name.endswith((".fastq", ".fq", ".gz")):
            file_path = os.path.join(target_folder, file_name)
            print(f" > Reading: {file_name}")
            
            counts, lens, ns = analyze_fastq(file_path)
            
            total_reads_count += len(lens)
            final_lengths.extend(lens)
            final_n_percents.extend(ns)
            
            for s, c in counts.items():
                final_counts[s] = final_counts.get(s, 0) + c
    
    # Sort and finish
    sorted_data = sorted(final_counts.items(), key=lambda x: x[1], reverse=True)
    save_top_sequences(sorted_data, report_name, total_reads_count)
    
    print(f"\nSUCCESS!")
    print(f"Total reads: {total_reads_count}")
    print(f"Report saved as: {report_name}")
    print("Opening graphs now...")
    
    create_plots(final_lengths, final_n_percents, target_folder)

else:
    print(f"\n[ERROR] Folder '{target_folder}' not found.")
    print("Make sure you typed it correctly. Available folders in this directory are:")
    # Show the user what folders are actually there to help them out
    for item in os.listdir():
        if os.path.isdir(item) and not item.startswith('.'):
            print(f" - {item}")