#!/usr/bin/python3
#-*- coding : utf-8 -*-

# AUTHORS = ("Hiba EL MOUSTANSIRI", "Alondra GONZÁLEZ")
# CONTACT = ("hiba.el-moustansiri@etu.umontpellier.fr","alondra.gonzalez-gonzalez@etu.umontpellier.fr")
# VERSION = "0.0.1"
# DATE = "16/12/2025"
# LICENCE = "This program is free software: you can redistribute it and/or modify
        #it under the terms of the GNU General Public License as published by
        #the Free Software Foundation, either version 3 of the License, or
        #(at your option) any later version.
        #This program is distributed in the hope that it will be useful,
        #but WITHOUT ANY WARRANTY; without even the implied warranty of
        #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        #GNU General Public License for more details.
        #You should have received a copy of the GNU General Public License
       # along with this program. If not, see <https://www.gnu.org/licenses/>."


# IMPORT MODULES #
import os, sys, re
import matplotlib.pyplot as plt

########################## 1 CHECK ##########################
def get_path(): # Function that asks the user for the path to a SAM file and verify some conditions to work with it.
    pathway = input("Enter the path to your SAM file: ")

    if not os.path.exists(pathway): # Check if the path exists
        print("Error: File doesn't exist") # If not: print error message
        return None # Stop function and return none in path.

    if not os.path.isfile(pathway): # Check if the path corresponds to a file
        print("Error: Path is not a file") # If not: print error message
        return None # Stop function and return none in path.

    if not pathway.lower().endswith('.sam'): # Check if file extension is .sam
        print("Error: File must have .sam extension") # If not: print error message
        return None # Stop function and return none in path.

    if os.path.getsize(pathway) == 0: # Check if file size is zero (empty file)
        print("Error: SAM file is empty") # If not: print error message
        return None  # Stop function and return none in path.

    return pathway # Return the validated file path
# MAIN PATH TEST


########################### 2 READ AND STOCKAGE ##########################
def mapq_threshold():
    answer = input("Do you want to apply a MAPQ threshold? (yes/no): ").strip().lower()

    if answer in ["yes", "y"]:
        while True:
            threshold = int(input("Enter the MAPQ threshold value: "))
            return threshold
    else:
        return None

def read_sam_file(path,mapq_threshold=None): # Define function to read SAM file
    dic_sam_data = {} # We create an empty dictionary to store reads
    line_read = 1 # Initialize our read counter starting by 1 (key for dictionary)

    file = open(path, 'r') # Open the SAM file in read mode

    # Read the first line in he sam file
    line = file.readline() # Readline: function that reads the first line of the sam file.

    while line:
        if not line.startswith("@"): # If the starts with @: we skip.
            columns = line.split() # Split the SAM line into columns by whitespace, because they're tab-separated.
            MAPQ = int(columns[4])
            if mapq_threshold is not None and MAPQ < threshold:  # Keep only reads above threshold
                continue
            # Variables that takes each position of the splited columns of sam file corresponding their information
            QNAME = columns[0]
            FLAG = columns[1]
            RDNAME = columns[2]
            POS = columns[3]
            CIGAR = columns[5]
            RNEXT = columns[6]
            PNEXT = columns[7]
            TLEN = columns[8]
            SEQ = columns[9]
            dic_sam_data[line_read] = { # Store all the read data
                "QNAME": QNAME,
                "FLAG": FLAG,
                "RNAME": RDNAME,
                "POS": POS,
                "MAPQ": MAPQ,
                "CIGAR": CIGAR,
                "RNEXT": RNEXT,
                "PNEXT": PNEXT,
                "TLEN": TLEN,
                "SEQ": SEQ
            }
            line_read += 1 # Increment read counter

        line = file.readline() # Move to next line of sam file

    file.close() # Close the SAM file
    return dic_sam_data # Return dictionary with all reads

# READING MAIN TEST



########################### 3 ANALYSE ##########################

#3.1 Convert the flag into binary ####
def flagBinary(FLAG) :
    flagBin = bin(int(FLAG)) # Transform the integer into a binary.
    flagBin = flagBin[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagBin = list(flagBin) # Convert into a list
    if len(flagBin) < 12: # Size adjustement to 12 (maximal flag size)
        add = 12 - len(flagBin) # We compute the difference between the maximal flag size (12) and the length of the binary flag.
        for i in range(add):
            flagBin.insert(0,'0') # We insert 0 to complete until the maximal flag size if it's < 12.
    return flagBin

# FLAG MAIN TEST
#print(flagBinary(5))

#3.2 Extract unmapped reads
# For quality control, we often want all unmapped reads in FASTA format.
# We also want a summary count for reporting all results.
def unmapped(dic_sam_data): # So, we create a fonction that identifies unmapped reads based on FLAG bits
    unmapped_count = 0 # Counter
    # We open two outputs:
    # - only_unmapped.fasta: stores sequences FASTA
    # - summary_file.txt: stores the final count of results
    # "a+" appends: to read and write, "w": to write.
    with open ("only_unmapped.fasta", "a+") as unmapped_fasta, open("summary_file.txt", "w") as summary_file:
        for key in dic_sam_data: # Iterate over all reads that we stored from the SAM file
            data_flag = flagBinary(dic_sam_data[key]["FLAG"]) # Convert numeric FLAG into binary list
            if int(data_flag[-3]) == 1:   # If the bit of the condition -3 is 1, the read is considered as "unmapped"
                unmapped_count += 1
                qname = dic_sam_data[key]["QNAME"] # QNAME becomes FASTA header to export into de unmapped fasta file
                seq = dic_sam_data[key]["SEQ"] # SEQ is the read's nucleotide sequence to export into de unmapped fasta file
                unmapped_fasta.write (f">{qname} function:unmapped\n{seq}\n") # Export the read in FASTA format
        summary_file.write("Total unmapped reads: " + str(unmapped_count) + "\n") # Write the unmapped count for the summary
    return unmapped_count

# 3.3 Extraire partially mapped reads:
# Some reads are mapped but with clipping/indels/skips (S/I/D/N/H etc.). These are not "fully matched".
# Here, "fully mapped" is defined very strictly as CIGAR == <number>M only.
# Anything else is considered "partially mapped".
def partiallyMapped(dic_sam_data):
    partially_mapped_count = 0 # Number of reads considered partially mapped
    mapped_count = 0 # Number of reads considered fully mapped (only M)
    with open ("only_partially_mapped.fasta", "a+") as partillay_mapped_fasta, open("summary_partially_mapped.txt", "w") as summary_file:
        for key, val in dic_sam_data.items():
            data_flag = flagBinary(val["FLAG"])
            if int(data_flag[-3]) == 0: # If the bit of the condition -3 read is 0 = mapped
                cigar = val["CIGAR"]
                # r"\d*[M]" Many valid full alignments are like "76M" (fully mapped), but soft-clipped reads "5S71M" become partial.
                if re.fullmatch(r"\d*[M]",cigar): # If CIGAR has only 'M', count as fully mapped
                    mapped_count +=1 # Increase fully mapped count
                else: # If there are any other, classify as partially mapped
                    partially_mapped_count += 1 # Increase partial count
                    partillay_mapped_fasta.write(f">{val['QNAME']}\n{val['SEQ']}\n") # Export sequence to FASTA
        summary_file.write("Total partially mapped reads: " + str(partially_mapped_count) + "\n"+"Mapped reads: "+ str(mapped_count)) # Save both counts
    return partially_mapped_count, mapped_count # Return both for later percentages

# 3.4 Run the complete complete analysis with:
# - counts unmapped, mapped and partially mapped reads
# - calculates their percentages
# - writes a summary file
# - returns all results so they can be reused outside the function

def full_mapping_analysis(dic_sam_data):
    unmapped_count = unmapped(dic_sam_data)
    partially_mapped_count, mapped_count = partiallyMapped(dic_sam_data)

    nb_total_reads = mapped_count + unmapped_count + partially_mapped_count # We compute the total number of reads analyzed

    mapped_percentage = (mapped_count / nb_total_reads) * 100 # Percentage of fully mapped reads
    unmapped_percentage = (unmapped_count / nb_total_reads) * 100  # Percentage of unmapped reads
    partially_mapped_percentage = (partially_mapped_count / nb_total_reads) * 100 # Percentage of partially mapped reads

    with open("summary_final_mapping_result.txt", "w") as summary_file: # Open the summary file in write mode
        summary_file.write(
            # Write percentages formatted with exactly two decimal places (.2f)
            f"Percentage of mapped reads: {mapped_percentage:.2f}%\n"
            f"Percentage of unmapped reads: {unmapped_percentage:.2f}%\n"
            f"Percentage of partially mapped reads: {partially_mapped_percentage:.2f}%\n"
            # Write the previous read counts
            f"Total number of reads: {nb_total_reads}\n"
            f"Mapped reads: {mapped_count}\n"
            f"Number of unmapped reads: {unmapped_count}\n"
            f"Number of partially mapped reads: {partially_mapped_count}\n"
        )
    # Return all results so they are accessible outside the function
    return (nb_total_reads, mapped_count, unmapped_count, partially_mapped_count,
            mapped_percentage, unmapped_percentage, partially_mapped_percentage)


# 3.5 Number of reads p/FLAG:
# 3.5.1 For paired-end sequencing, a useful QC is to count cases where one read maps but its mate does not.
def Pair_mapped_non_mapped(dic_sam_data):
    #mapq_min = input("Entrez le score de mapping : ")
    pair_mapped_unmapped_count = 0 # Counter for the number of pairs
    with open ("summary_pairs.txt", "w") as summary_file: # Save results to file
        for key, val in dic_sam_data.items():
            data_flag = flagBinary(val["FLAG"])  # Convert FLAG to binary for paired-end bit tests
            mapq = int(val["MAPQ"])  # Convert MAPQ to int.
            if int(data_flag[-1]) == 1 and int(data_flag[-4]) == 1: # Verify that the read are "paired" and "mate unmapped"
                cigar = val["CIGAR"] # We get CIGAR to ensure the read itself is fully mapped
                if re.fullmatch(r"\d*[M]",cigar): # If read is fully mapped (only M)
                    pair_mapped_unmapped_count +=1 # Increase counter
        summary_file.write("Total pairs with one well-mapped read and the other unmapped: "+str(pair_mapped_unmapped_count))
    return pair_mapped_unmapped_count

# 3.5.2 Calculate the Number of read pairs with one fully mapped read and the other partially mapped
def Pair_mapped_partiel_mapped(dic_sam_data):
    pair_mapped_partiel_mapped_count = 0 # Counter
    with open ("summary_pairs.txt", "w") as summary_file:
        for key, val in dic_sam_data.items():
            data_flag = flagBinary(val["FLAG"]) # Convert FLAG for bit checking
            mapq = val["MAPQ"]  # Convert MAPQ to int
            if int(data_flag[-2]) == 1 and int(data_flag[-3]) == 0: # Verify that te read is paired and mapped
                cigar = val["CIGAR"]
                if not re.fullmatch(r"\d*[M]",cigar): # If it is not fully mapped, treat as partially mapped
                    pair_mapped_partiel_mapped_count +=1 # Increase count
        summary_file.write("Total des paires avec un read totalement mappé et l'autre partiellement mappé: "+str(pair_mapped_partiel_mapped_count)) # Save output
    return pair_mapped_partiel_mapped_count

# 4. With what quality are reads mapped?
# MAPQ is a confidence score. Counting reads by MAPQ helps you see if most alignments are high-confidence or not.
def qualityCount(dic_sam_data):
    quality_reads = {}
    for key, val in dic_sam_data.items(): # Iterate over reads
        data_flag = flagBinary(val["FLAG"]) # Convert FLAG
        if int(data_flag[-3]) == 1:  # If read is unmapped, MAPQ is not meaningful for alignment confidence so we skip
            continue # Skip to next read
        mapq = val["MAPQ"]  # Convert quality to integer.
        if mapq not in quality_reads: # If first time seeing this MAPQ score, initialize
                quality_reads[mapq] = 1 # Start count at 1
        else:
                quality_reads[mapq] += 1 # Otherwise increment count
    return quality_reads


# 3.6 Chromosome position analysis
# We often want to know how reads distribute across chromosomes, and where on each chromosome they map.
def analyze_chromosome_positions(dic_sam_data):
    chromosome_positions = {}  # Dictionary to store positions per chromosome.
    for key, val in dic_sam_data.items(): # Iterate all reads
        data_flag = flagBinary(val["FLAG"]) # Convert FLAG to binary
        if int(data_flag[-3]) == 0: #Only process mapped reads
                chrom = val["RNAME"] # Chromosome/contig name
                pos = int(val["POS"]) # Convert position to int so it can be sorted/aggregated numerically
                if chrom not in chromosome_positions:  # If first time we see this chromosome, create an empty list
                    chromosome_positions[chrom] = [] # Initialize list for that chromosome
                chromosome_positions[chrom].append(pos)  # Add the read's mapping position to that chromosome list
    return chromosome_positions


############################# MAIN ###########################################
path = get_path() # Call the function and store the valid path
print("Valid path")

threshold=mapq_threshold() #store the mapq threshold
dic_sam_data = read_sam_file(path) # Read and store SAM data

### Print the summary statistics for mapped reads
print("Nombre de reads non mappés :", unmapped(dic_sam_data))
print("Number of partially mapped reads :", partiallyMapped(dic_sam_data))
(nb_total_reads, mapped_count, unmapped_count, partially_mapped_count,
 mapped_percentage, unmapped_percentage, partially_mapped_percentage) = full_mapping_analysis(dic_sam_data)
print("Final mapping results written to summary_final_mapping_result.txt")
print("Number of mapped/unmapped read pairs:", Pair_mapped_non_mapped(dic_sam_data))
print("Nombre de paires de reads partiellement mappé/-mappés :", Pair_mapped_partiel_mapped(dic_sam_data))
print("Nombre des reads selon leur qualité:", qualityCount(dic_sam_data))

# Bar graph to represent the quality distribution (MAPQ).
quality_reads=qualityCount(dic_sam_data)
x_quality = []  # List of quality values (MAPQ).
y_quality = []  # List for the number of reads for each quality.

for x, y in quality_reads.items():   # Scans the quality/number pairs.
    x_quality.append(x)  # Adds the quality to the X list.
    y_quality.append(y)   # Adds the number of reads to the Y list.

plt.bar(x_quality, y_quality, color='skyblue', edgecolor='black')  # Creates the bar graph.
plt.xlabel('Mapping quality')  # Adds a label for the X axis.
plt.ylabel('Number of reads')  # Adds a label for the Y axis.
plt.yscale("log") # Sets the scale of the Y axis to logarithmic to better visualise large differences.
plt.title(f'Number of reads by mapping quality')  # Adds a title.
plt.savefig("mapping_quality.png", dpi=600, bbox_inches="tight") # Saves the plot.
plt.clf()  # Cleans up the figure.

chromosome_positions=analyze_chromosome_positions(dic_sam_data)
#print("Nombre des reads par chromosome:", analyze_chromosome_positions(dic_sam_data))