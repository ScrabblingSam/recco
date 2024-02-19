#!/usr/bin/env python
# coding: utf-8

from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import tempfile
import os
import shutil
import argparse
import sys
import glob
#take arguments from terminal
def take_args():
    parser=argparse.ArgumentParser(description="ReCCO.py Version: 1.0, Developer: Peijie Yan, Kyoto University. -Rearrange and Concatenate fasta files in order(blast+ and Biopython installation required and please make sure the fasta files end with .fasta, .fa, or .fna)")
    
    parser.add_argument("-i", "--input_dir", required=True, help="Required: Directory of input file containing both query and reference fasta files")
    parser.add_argument("-r", "--reference_file_name", required=True, help="Required: Name of your reference file")
    parser.add_argument("-n", "--number_of_N", required=False, type=int, default=1000, help="Number of Ns between concatenated contigs. Default value: 1000")
    parser.add_argument("-o", "--output_dir", required=True, help="Required: Directory for output files")
    parser.add_argument("-s", "--suffix", default="_ReCCO", help="Suffix for rearranged sequences. Default: '_ReCCO'")
    parser.add_argument("-m", "--min_length", required=False, type=int, default=0, help="Minimum length of input contigs retained in the output. Default value: 0")
    parser.add_argument("-b", "--blast_bin_dir", required=False, help="Directory for blast bin folder")
    
    try:
        args=parser.parse_args()
        query_dir=os.path.abspath(args.input_dir)
        reference_name=os.path.basename(args.reference_file_name)
        work_space=os.path.abspath(args.output_dir)
        mini_length=args.min_length
        number_of_N=args.number_of_N
        suffix=args.suffix
        
        if not args.blast_bin_dir:
            blastn=subprocess.run("which blastn", shell=True, text=True, capture_output=True)
            if blastn.returncode ==0:
                blastn_path=blastn.stdout.strip()
                blast_bin_dir=os.path.dirname(blastn_path)
        
            else:
                raise ValueError("blast_bin_dir argument not provided and failed to find blast bin dir in environment variable.")
        else:
            blast_bin_dir = os.path.abspath(args.blast_bin_dir)
        
        
        
        if os.path.exists(work_space):
            while True: 
                choice = input("Warning: The output directory already exists. Do you wish to backup the existing directory before execution? <y/n> or type q to exit:").lower()
                if choice == "y":
                #backup the data and clear the output directory
                    try:
                        backup_dir_base = input("Please enter the directory at which you want to back up the data:")
                        counter=1
                        backup_dir = backup_dir_base+"/backup_data_ReCCO"
                        while os.path.exists(backup_dir):
                            counter+=1
                            backup_dir=f"{backup_dir_base}/backup_data_ReCCO{counter}"
                        shutil.copytree(work_space, backup_dir)
                        print("The backup data is created in", backup_dir)
                        for filename in os.listdir(work_space):
                            file_full_path=os.path.join(work_space,filename)
                            if filename!="backup_data":
                                if os.path.isfile(file_full_path) or os.path.islink(file_full_path):
                                    os.unlink(file_full_path)
                                elif os.path.isdir(file_full_path):
                                    shutil.rmtree(file_full_path)
                        break  
                    except Exception as e:
                        print(f"Error while backing up data: {e}")
                        sys.exit(1)
                elif choice == "n":
                    choice2 = input("If the program continues, all files in the given output file directory will be deleted. This action cannot be undone. Do you want to continue? <y/n>").lower()
                    if choice2 == "y":
                        try:
                            shutil.rmtree(work_space)
                            os.makedirs(work_space)
                            break 
                        except OSError as e:
                            print(f"Error while cleaning output directory: {e}")
                            sys.exit(1)
                    elif choice2 == "n":
                        sys.exit(0)
                    else:
                        print("Invalid input. Please enter y or n.")  
                elif choice == "q":
                    sys.exit(0)
                else:
                    print("Invalid input. Please enter y, n, or q to quit.")  
        else:
            os.makedirs(work_space)

        
        small_chain_files=[]
        input_file_main_chain=None
        
        for filename in os.listdir(query_dir):
            if filename.endswith(".fasta") or filename.endswith(".fa") or filename.endswith(".fna") :
                file_path=os.path.join(query_dir, filename)
                if filename==reference_name:
                    input_file_main_chain=file_path
                else:
                    small_chain_files.append(file_path)
        
    except Exception as e:
        if e.code==0:
            sys.exit(0)
        else:
            print("Error: Argument parsing failed. Please check if there are spaces between query fasta file directories")
            sys.exit(1)
    
    return small_chain_files, input_file_main_chain, blast_bin_dir, work_space, mini_length, number_of_N, suffix, reference_name

small_chain_files, input_file_main_chain, blast_bin_dir, work_space, mini_length, number_of_N, suffix, reference_name=take_args()

#Read the reference fasta file and if the reference file contains multiple sequences, concatenate them into one contig
def read_fasta_file_main(input_file_main_chain):
    main_chain_sequences=[]
    for record in SeqIO.parse(input_file_main_chain,"fasta"):
        main_chain_sequences.append(str(record.seq))
#create a new file so that you don't erase the original data
    reference_path=os.path.join(work_space, "reference_contig.fasta")
    with open(reference_path, 'w') as file:
        file.write(">"+os.path.basename(input_file_main_chain).split('.')[0]+"_CR"+'\n')
        file.write("".join(main_chain_sequences))
    return reference_path
reference_path=read_fasta_file_main(input_file_main_chain)

#Read the query chain fasta files and count the length of each seqeunce and throw away very small contigs
def read_fasta_file_small(small_chain_files, mini_length):
    sequences={}
    for file in small_chain_files:
        filename=os.path.basename(file.rstrip('/'))
        basename, extension=os.path.splitext(filename)
        filename=basename+".fasta"
        file_sequences={}
        for record in SeqIO.parse(file,"fasta"):
            sequence_id=record.id
            sequence=str(record.seq)
            if len(sequence)>mini_length:
                file_sequences[sequence_id]=sequence
        sequences[filename]=file_sequences
    return sequences
sequences=read_fasta_file_small(small_chain_files, mini_length)

#This is to run the blast search. 
#database made of reference chain
def runblast(input_file_main_chain, sequences):
    db_name = "reference_database"
    db_path= os.path.join(work_space, db_name)
    makeblastdb_cmd = f"{blast_bin_dir}/makeblastdb -in {reference_path} -dbtype nucl -out {db_path}"
    subprocess.run(makeblastdb_cmd, shell = True, check = True)
#run blast against the reference chain and get blast results files, the number of which is the same as the number of query contigs
    blast_results={}
    for filename in sequences:
        blast_result=[]
        for sequence_id, sequence in sequences[filename].items():
            print("Running Blast for sequence:", sequence_id)
            with tempfile.NamedTemporaryFile(mode='w') as fp:
                fp.write(f'>{sequence_id}\n{sequence}\n')
                fp.seek(0)
                blast_cmd = f"{blast_bin_dir}/blastn -query {fp.name} -db {db_path} -outfmt '6 std qlen slen'"
                result = subprocess.run(blast_cmd, shell = True, check = True, capture_output = True, text = True)
                blast_result.append(result.stdout)
        blast_results[filename]=blast_result
    print("Blast search completes")
    return blast_results
blast_results=runblast(input_file_main_chain, sequences)

try:
    db_files_pattern = os.path.join(work_space, "reference_database.n*")

    db_files_list = glob.glob(db_files_pattern)
    for db_file in db_files_list:
        os.remove(db_file)
except OSError as e:
    print(f"Error while cleaning main_contig_database: {e}")
    
#write the blast results out to files in a folder called blast_results
blast_results_dir=os.path.join(work_space, "blast_results")
if os.path.exists(blast_results_dir):
    shutil.rmtree(blast_results_dir)
os.makedirs(blast_results_dir, exist_ok=True)
for i, filename in enumerate(blast_results):
    filename_no_extension=os.path.splitext(filename)[0]
    out_put_blast_path=os.path.join(blast_results_dir, filename_no_extension)
    with open(out_put_blast_path, "w") as output_handle:
        output_handle.write("".join(blast_results[filename]))

#To filter the outfmt 6 files and to reverse compliment sequences if any
def edit_outfmt6_file(blast_results_dir):
    try:
        for blast_result_file in os.listdir(blast_results_dir):
            sequences_key=blast_result_file+".fasta"
            if blast_result_file=='.DS_Store':
                continue
            blast_result_file_path=os.path.join(blast_results_dir, blast_result_file)
            #make sure the dirctory points to a file
            if os.path.isfile(blast_result_file_path):
                best_alignment={}
                outfmt6_results=[]
                reversed_sequence_id=[]
                with open(blast_result_file_path, "r") as file:
                    for line in file:
                        datas=line.strip().split('\t')
                    #Update dictionary with longest alignment length for each query
                        query_id=datas[0]
                        subject_id=datas[1]
                        align_length=int(datas[3])
                        e_value=float(datas[10])
                        pident=float(datas[2])
                        s_start=int(datas[8])
                        s_end=int(datas[9])
                        if query_id in best_alignment:
                            if align_length>best_alignment[query_id][-1][1]:
                                best_alignment[query_id]=[]
                                best_alignment[query_id].append([subject_id, align_length, e_value])
                                #remove last 14 datas elements
                                for i in range(14):
                                    outfmt6_results.pop()
                            else:
                                continue
                        else:
                            if e_value<0.1:
                                best_alignment[query_id]=[]
                                best_alignment[query_id].append([subject_id, align_length, e_value, pident])
                            else:
                                continue
                        #if the order of the segment is reverse complemented
                        if s_start>s_end:
                            if query_id in sequences[sequences_key] and query_id not in reversed_sequence_id:
                                sequence=sequences[sequences_key][query_id]
                                sequence=Seq(sequence).reverse_complement()
                                sequences[sequences_key][query_id]=sequence
                                reversed_sequence_id.append(query_id)
                                e=datas[8]
                                datas[8]=datas[9]
                                datas[9]=e


                        outfmt6_results.append(datas)
                    #Sort blast results based on subject start position
                    sorted_outfmt6_results = sorted(outfmt6_results, key=lambda x: int(x[8]))
                with open(blast_result_file_path, "w") as file:
                    for result in sorted_outfmt6_results:
                        file.write('\t'.join(result) + '\n')
    except OSError:
        print("Error in opening blast results file. Please check blast is run properly.")
        sys.exit(1)

edit_outfmt6_file(blast_results_dir)
     
# To edit the fasta file
# Based on filtered blast file, concatenate the sequences.
def edit_fasta_file(sequences):
    #put sequences that are in sequences dictionary but not in outfmt6 file to the end. 
    edited_sequence_id=[]
    query_subject={}
    subjects={}
    for filename in os.listdir(blast_results_dir):
        sequences_to_end={}
        sequences_to_end_id=[]
        if filename == '.DS_Store':#to skip mac generated system file
            continue
        filename_fasta=os.path.splitext(filename)[0]+".fasta"
        blast_result_file_path=os.path.join(blast_results_dir, filename)
        with open(blast_result_file_path, "r") as file:
            #To capture cases when the edited blast results files are empty. 
            content=file.readline()
            if not content:
                print(f"The blast result file {filename} is empty. Skipping...")
                continue
            file.seek(0)
            for sequence_id in sequences[filename_fasta]:
                found_sequence = False  # Flag to track if sequence is found in blast results
                for line in file:
                    datas = line.strip().split('\t')
                    query_id = datas[0]
                    if sequence_id == query_id:
                        found_sequence = True
                        break
                file.seek(0)
                if not found_sequence:
                    sequences_to_end_id.append(sequence_id)
            #don't delete the sequences. Put them to the end of the contig
            for sequence_id in sequences_to_end_id:
                sequences_to_end[sequence_id]=sequences[filename_fasta][sequence_id]
                del sequences[filename_fasta][sequence_id]


        with open(blast_result_file_path, "r") as output_handle:
            for line in output_handle:
                fields=line.strip().split('\t')
                #percentage of identical matches
                pident=float(fields[2])
                query_id=fields[0]
                subject_id=fields[1]
                align_length=int(fields[3])
                s_start=int(fields[8])
                s_end=int(fields[9])
                if query_id not in query_subject:
                    query_subject[query_id]=[]
                query_subject[query_id].append([query_id, subject_id, s_start, s_end, align_length])
                #reorder the small sequences  
                if subject_id in subjects:
                    query_id_before=subjects[subject_id][-1]
                    if query_id in sequences[filename_fasta]: #and query_id not in edited_sequence_id:
                        sequence=sequences[filename_fasta][query_id]
                        #situation 1 
                        if query_subject[query_id][0][2]>=query_subject[query_id_before][0][3]:
                            sequence=sequences[filename_fasta][query_id_before]+"N"*number_of_N+sequence
                            sequences[filename_fasta][query_id]=sequence
                            edited_sequence_id.append(query_id_before)
                            edited_sequence_id.append(query_id)
                            subjects[subject_id].append(query_id)
                            del sequences[filename_fasta][query_id_before]
                        #also handles situations where the start_positions are the same
                        elif query_subject[query_id][0][2]<query_subject[query_id_before][0][3] or query_subject[query_id_before][0][2]==query_subject[query_id][0][2]:
                            # c_id=(query_subject[query_id][0][2]+query_subject[query_id][0][3])/2
                            # c_id_before=(query_subject[query_id_before][0][2]+query_subject[query_id_before][0][3])/2
                            # commented out section tries to sort the query contigs based on the central points.
                        #situation 2: Sort by start positions 
                            sequence=sequences[filename_fasta][query_id_before]+sequence
                            sequences[filename_fasta][query_id]=sequence
                            edited_sequence_id.append(query_id_before)
                            edited_sequence_id.append(query_id)
                            subjects[subject_id].append(query_id)
                            del sequences[filename_fasta][query_id_before]
                            # #situation 3 
                            # elif c_id < c_id_before:
                            #     sequence=sequence+sequences[filename_fasta][query_id_before]
                            #     sequences[filename_fasta][query_id]=sequence
                            #     edited_sequence_id.append(query_id_before)
                            #     edited_sequence_id.append(query_id)
                            #     subjects[subject_id].append(query_id)
                            #     del sequences[filename_fasta][query_id_before]
                        
                    else:
                        print("error in editing fasta file")
                else :
                        subjects[subject_id]=[]
                        subjects[subject_id].append(query_id)            
                                              
        keys=list(sequences[filename_fasta].keys())
        for i in range(len(keys)):
            sequences[filename_fasta][filename+suffix]=sequences[filename_fasta][keys[i]]
            del sequences[filename_fasta][keys[i]]
        for ids in sequences_to_end:
            sequences[filename_fasta][filename+suffix]+=("N"*number_of_N+sequences_to_end[ids])
            sequences[filename_fasta][filename+suffix]="".join(sequences[filename_fasta][filename+suffix])  
        keys2=list(subjects.keys())
        for i in range(len(keys2)):
            subjects[filename+suffix]=subjects[keys2[i]]
            del subjects[keys2[i]]
    return sequences,edited_sequence_id,subjects

sequences,concatenated_sequence_id,subjects=edit_fasta_file(sequences)  

def write_to_new_fasta_file(sequences):
    edited_fasta_files_dir=os.path.join(work_space, "edited_fasta_file")
    if os.path.exists(edited_fasta_files_dir):
        shutil.rmtree(edited_fasta_files_dir)
    os.makedirs(edited_fasta_files_dir, exist_ok=True)
    for filename in sequences:
        filename_no_extension=os.path.splitext(filename)[0]
        output_fasta_path=os.path.join(edited_fasta_files_dir, filename_no_extension+suffix+".fasta")
        with open(output_fasta_path,"w") as new_fasta:
            for sequence_id,sequence in sequences[filename].items():
                new_fasta.write(f'>{sequence_id}\n')
                new_fasta.write(f'{sequence}\n')
    return edited_fasta_files_dir

edited_fasta_files_dir=write_to_new_fasta_file(sequences)

#concatenate the results and reference chain into one fasta file so that it can be passed to digalign
def give_results():
    results_path=os.path.join(work_space,"ToDiGAlign.fasta")
    with open(results_path, "w") as output_handle:
        for filename in os.listdir(edited_fasta_files_dir):
            if filename.endswith(".fasta"):
                fasta_sequences=SeqIO.parse(os.path.join(edited_fasta_files_dir,filename), 'fasta')
                for record in fasta_sequences:
                    SeqIO.write(record, output_handle, "fasta")
        reference_sequences=SeqIO.parse(reference_path,"fasta")
        for record2 in reference_sequences:
            SeqIO.write(record2, output_handle, "fasta")
give_results()
#rewrite the format of blast_results to tsv and rename the reference file. 
def increase_usability(blast_results_dir):
    shutil.move(reference_path, os.path.join(edited_fasta_files_dir, os.path.splitext(reference_name)[0]+"_CR.fasta"))
    for filename in os.listdir(blast_results_dir):
        if not filename.endswith(".tsv"):
            new_filename=filename+".tsv"
            os.rename(os.path.join(blast_results_dir,filename),os.path.join(blast_results_dir,new_filename))
increase_usability(blast_results_dir)
print("ReCCO completes")