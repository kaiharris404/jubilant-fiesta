import os
import gzip
from Bio import SeqIO
import matplotlib.pyplot as plt

# ==========================================
# MODULE 1: Tool 1 & 2 (Filtering and Counting)
# ==========================================
def analyze_fastq(filename):
    """
    Reads the file, throws out bad reads (>5% N), 
    and collapses duplicates into a dictionary.
    """
    seq_counts = {}
    lengths = []
    n_percents = []

    # Open compressed or uncompressed files
    if filename.endswith(".gz"):
        handle = gzip.open(filename, "rt")
    else:
        handle = open(filename, "r")

    for record in SeqIO.parse(handle, "fastq"):
        seq = str(record.seq).upper()
        seq_len = len(seq)
        
        # Prevent math errors if a sequence is completely empty
        if seq_len == 0:
            continue
            
        # Calculate N% = (N count / length) * 100
        n_count = seq.count('N')
        n_percent = (n_count / seq_len) * 100
        
        # --- TOOL 2: N-Content Filtering ---
        # If it has more than 5% Ns, we skip it entirely!
        if n_percent > 5.0:
            continue 
            
        # --- TOOL 1: Redundancy Reduction ---
        # If it survived the filter, add it to our dictionary
        seq_counts[seq] = seq_counts.get(seq, 0) + 1
            
        # Keep track of metrics for our charts
        lengths.append(seq_len)
        n_percents.append(n_percent)

    handle.close() 
    return seq_counts, lengths, n_percents

# ==========================================
# MODULE 2: Ranking
# ==========================================
def sort_by_abundance(seq_dict):
    """Ranks sequences by abundance (highest count first)."""
    # Sorting a dictionary by its values
    return sorted(seq_dict.items(), key=lambda x: x[1], reverse=True)

# ==========================================
# MODULE 3: Outputs (Dashboard, FASTA, Charts)
# ==========================================
def save_reports(sorted_seqs, text_file, fasta_file, total_reads):
    """Saves the abundance table (.txt) and deduplicated FASTA (.fasta)."""
    
    # 1. Save the Quality Dashboard / Abundance Table
    with open(text_file, 'w') as t_file:
        t_file.write("--- QUALITY DASHBOARD & ABUNDANCE TABLE ---\n")
        t_file.write(f"Total reads passing 5% N filter: {total_reads}\n")
        t_file.write(f"Unique sequences found: {len(sorted_seqs)}\n\n")
        
        for i, (seq, count) in enumerate(sorted_seqs[:50]): # Top 50 for the table
            t_file.write(f"Rank {i+1} | Count: {count} | Seq: {seq}\n")

    # 2. Save the Deduplicated FASTA
    with open(fasta_file, 'w') as f_file:
        for i, (seq, count) in enumerate(sorted_seqs):
            # FASTA format requires a > header line before the sequence
            f_file.write(f">unique_seq_{i+1}_abundance_{count}\n")
            f_file.write(f"{seq}\n")

def create_plots(lengths, n_percents, sorted_seqs, folder_name):
    """Generates the bar chart and histograms requested in the rubric."""
    
    # Make a wide canvas for 3 charts
    plt.figure(figsize=(15, 5))
    
    # CHART 1: Abundance Bar Chart (Top 10)
    plt.subplot(1, 3, 1)
    top_10 = sorted_seqs[:10]
    # Create labels like "Seq 1", "Seq 2" so they fit on the screen
    x_labels = [f"Seq {i+1}" for i in range(len(top_10))]
    counts = [item[1] for item in top_10]
    
    plt.bar(x_labels, counts, color='lightgreen', edgecolor='black')
    plt.title('Top 10 Sequence Abundance')
    plt.xticks(rotation=45) # Tilt labels so they don't overlap
    plt.ylabel('Read Count')
    
    # CHART 2: Read Length Histogram
    plt.subplot(1, 3, 2)
    plt.hist(lengths, bins=30, color='skyblue', edgecolor='black')
    plt.title('Read Length Distribution')
    plt.xlabel('Length (bp)')
    plt.ylabel('Frequency')
    
    # CHART 3: N-Content Histogram
    plt.subplot(1, 3, 3)
    plt.hist(n_percents, bins=30, color='salmon', edgecolor='black')
    plt.title('N-Content Distribution (<5%)')
    plt.xlabel('% of N bases')
    
    plt.tight_layout()
    plt.show()

# ==========================================
# MAIN EXECUTION BLOCK (Interactive)
# ==========================================
print("--- Preprocessing Tools (Redundancy & N-Filtering) ---")
target_folder = input("Enter the folder name to analyze (e.g., barcode09): ").strip()

if os.path.exists(target_folder) and os.path.isdir(target_folder):
    
    # Setup our two output file names
    dashboard_file = f"{target_folder}_dashboard.txt"
    fasta_out = f"{target_folder}_deduplicated.fasta"
    
    final_counts = {}
    final_lengths = []
    final_n_percents = []
    total_good_reads = 0

    print(f"\nProcessing files inside '{target_folder}'...")

    for file_name in os.listdir(target_folder):
        if file_name.endswith((".fastq", ".fq", ".gz")):
            file_path = os.path.join(target_folder, file_name)
            print(f" > Reading and filtering: {file_name}")
            
            counts, lens, ns = analyze_fastq(file_path)
            
            total_good_reads += len(lens)
            final_lengths.extend(lens)
            final_n_percents.extend(ns)
            
            for s, c in counts.items():
                final_counts[s] = final_counts.get(s, 0) + c
    
    # Run the final ranking and saving
    sorted_data = sort_by_abundance(final_counts)
    save_reports(sorted_data, dashboard_file, fasta_out, total_good_reads)
    
    print(f"\nComplete!")
    print(f"Total reads passing filter: {total_good_reads}")
    print(f"Quality Dashboard saved: {dashboard_file}")
    print(f"Deduplicated FASTA saved: {fasta_out}")
    print("Opening the 3 requested charts now...")
    
    create_plots(final_lengths, final_n_percents, sorted_data, target_folder)

else:
    print(f"\n[ERROR] Folder '{target_folder}' not found.")
