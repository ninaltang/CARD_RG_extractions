"""
This script will extract resistance gene sequences in FASTA format.
You must have the CARD output for this isolate in csv format and the full
genome FASTA sequence for each isolate in the same folder

Part 1 (lines 8-37): Set path to find relevant files
Part 2 (lines 38-63): Extracting resistance gene sequence
Part 3 (lines 65-90): Writing forward sequence of resistance gene
"""
import os
samples_li = ["NR6283"] #enter samples in list for faster computing
for i in samples_li:
    sample = i
    og_path = r"sample_path/{}".format(sample) #enter path for sample folder
    os.chdir(og_path)
    sample_card = "{}_CARD.csv".format(sample) #open CARD csv
    card_csv = open(sample_card,"r")
    line1 = card_csv.readline() #read heading line
    line1 = card_csv.readline() #read actual content line
    while line1 != "":
        line1 = line1.split(",") 
        line1= tuple(line1)
        isolate_name = str(line1[1][:6])
        contig_name = str(line1[1])
        genome_name = contig_name[0:6]
        rg_name = str(line1[8])
        rg_start = int(line1[2])-1
        rg_end = int(line1[3])-1
        temp_contig_name = list(contig_name)
        final_contig_name = str(contig_name[0:7])
        end_contig_name = str(contig_name[7:])
        for c in end_contig_name:
            if c != "_":
                final_contig_name = final_contig_name + c  
            else:
                break
        final_contig_name = "".join(final_contig_name)
        fas_contig_name = ">" + final_contig_name #fas_contig_name is the final_contig_name with > in front which is fasta format
        genome_name = "{}.fasta".format(genome_name) #name of genome fasta file
        genome_file = open(genome_name,"r")
        final_genome = ""
        genome_line1 = genome_file.readline() #reads first line of genome_line1 which is the header
        #print("genome header is ",genome_line1) #use this as a check if having problems
        while genome_line1 != "": #when this line has some content, we enter this loop    
            genome_line1 = genome_line1.split()
            #print("splitting header") #use this as a check if having problems
            if genome_line1[0] == fas_contig_name: #if the first object in the list genome_line1 which was the header separated into a list at the spaces, is the same as the fas_contig_name, then enter loop
                genome_line1 = genome_file.readline() #read the next line which is the start of the data
                #print("reading first line of correct contig")
                while genome_line1 != "":
                    if genome_line1[0] == ">": #while the first index of this next line (string) is not another header, enter this loop
                        #print("finished reading contig")
                        break
                    else:
                        genome_line1 = genome_line1.strip() #get rid of the \n at the end
                        #print("reading contig")
                        final_genome = final_genome + genome_line1 #append this to the final_genome list
                        #print("adding line to final genome")
                        genome_line1 = genome_file.readline() #read the next line
                        #print("reading next line")
            else: #if the first object in the list genome_line1 which was the header separated into a list at the spaces is not the same as the fas_contig_name, then enter this loop
                genome_line1 = genome_file.readline() #keep reading lines until we either hit the contig name or we run out of lines of the file
        #print("final genome finished")
        rg_dna = final_genome[rg_start:(rg_end+1)]
        #print("final resistance gene sequence: ",rg_dna)
        rg_ori = str(line1[4])
        new_path = r"{}/gene_extractions".format(og_path)
        try:
            os.mkdir(new_path)
        except OSError:
            print("Creation of the directory {} failed".format(new_path))
        else:
            print("Successfully created the directory {}.".format(new_path))
        os.chdir(new_path)
        rg_file = open("{}_{}".format(contig_name,rg_name), "a")
        if line1[4] == "+":
            new_line = ">{} {} \n"
            new_line = new_line.format(contig_name,rg_name)
            rg_file.write(new_line)
            rg_file.write(rg_dna)
        elif line1[4] == "-":
            new_line = ">{} {} (reverse) \n"
            new_line = new_line.format(contig_name,rg_name)
            rg_file.write(new_line)
            rg_file.write(rg_dna)
        rg_file.close()
        line1 = card_csv.readline()
        os.chdir(og_path)
        #print("Finished gene extraction. Moving onto next gene")
    card_csv.close()
    print("Finished all extractions.")
