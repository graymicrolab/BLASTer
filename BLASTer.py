# BLASTer
version = "4"
updated = "April 20, 2020"

"""
BLASTer is a tool which automates a workflow for identifying homologs of
proteins, examining their conservation across different species, and
generating useful alignments and tables for publication. It is basically
automation of what I want to know about any protein I plan on working with.

Based on Python 2 code previously found in scripts BLASTer_v3.3.py and
align-o-tron2002.py. (Both of which date from ~2016.)

Python 3.8.2
Biopython 1.76
matplotlib 3.2.1
MUSCLE 3.8.31 (www.drive5.com/muscle/downloads.htm)
BLAST+ 2.10.0 (ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
Logomaker 0.8 (Tareen & Kinney 2019 Bioinformatics btz921. bioRxiv doi:10.1101/635029)

Setup (saved to file; adjustable if needed):
    Inputs:
        Entrez email
        save path
        path to standalone MUSCLE and blastp applications

Step 1: Initial Inputs
    Inputs:
        protein name
        protein sequence (full or partial)

Step 2: The Align-o-Tron
    Outputs:
        (optional) graph of hits vs. e-value
        (optional) automatically generated alignments and sequence logos at
            e-values where this makes sense
    
Step 3: BLAST Parameters (if you choose to do higher quality alignment)
    Inputs:
        BLAST e-value threshold
        maximum number of hits to return
    Outputs:
        temporary XML file of raw BLAST results
        number of hits
        tables of species and genus occurence
        (optional) repeat at a different e-value threshold

Step 4: Pruning the List & Converting the Output File
    Inputs:
        pruning?
        if so, by species or genus?
        how many per species / genus?
    Outputs:
        FASTA-formatted list of full sequences of homologous proteins (delete XML file)
        number of proteins in that list

Step 5: Trim Contents of Protein List By Length?
    Outputs:
        length distribution (relative to the full length of the top hit)
    Inputs:
        trimming?
        if so, how and by how much?
    Outputs:
        trimmed FASTA-formatted protein list (overwrites untrimmed file)
        number of proteins in that list

Step 6: Generate Alignment
    Outputs: (optional)
        MUSCLE-aligned, FASTA-formatted alignment file
        sequence logo

Step 7: Generate Homology Table
    Outputs: (optional) 
        table of homology and taxonomy
"""

# MODULES

# standard Python modules
import os
import re
import sys
import math
import datetime

# Biopython modules
from Bio import SeqIO
from Bio import Entrez
from Bio import AlignIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio.Blast.Applications import NcbiblastpCommandline

# modules with specialized functions
import numpy as np
import pandas as pd
import logomaker as lm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# HELP SYSTEM TOPICS

help_dict = {
    "license": "Copyright 2020 by Michael J. Gray. This work is licensed under the GNU GPLv3 License. To view a copy of this license, visit choosealicense.com/licenses/gpl-3.0/",
    "contact": "You can contact the author of this software at mjgray@uab.edu.",
    "tech support": "Unfortunately, due to my time constraints, this software is available ‘as-is’ and I cannot promise any tech support or future updates. I hope it’s helpful, but it may need some tweaking, especially as the underlying tools get updated in the future.",
    "credits": "BLASTer was written with Python 3.8.2, Biopython 1.76, MUSCLE 3.8.31, BLAST+ 2.10.0, matplotlib 3.2.1, and Logomaker 0.8.",
    "dependencies": "BLASTer requires Python 3 and the following Python modules: biopython 1.76 or greater, matplotlib 3.2.1 or greater, and logomaker 0.8 or greater, all available for installation with pip.",
    "tutorial": "A tutorial is available for download at https://github.com/graymicrolab/BLASTer",
    "Entrez email": "NCBI requires a valid email address to make automated requests from their Entrez server.",
    "standalones": "BLASTer requires the MUSCLE 3.8.31 executable (available from www.drive5.com/muscle/downloads.htm) and the blastp exectutable from BLAST+ 2.10.0 (available from ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) in the folder defined in the settings file.",
    "save path": "BLASTer will create folders containing your results in the folder defined in the settings file.",
    "protein name": "BLASTer will accept any text as the name of your protein of interest, and that name will be used to generate the names of the resulting output files, along with a datestamp in the format YEARMODY (e.g 20200408 for April 8, 2020).",
    "overwriting data": "If a folder with the name of your requested analysis already exists, BLASTer can overwrite it. It will give you the option to quit and change the name of that folder instead, though.",
    "sequence": "You can enter either a partial or a full-length protein sequence. BLASTer will use whatever sequence you enter for its BLAST search, but will return full-length protein sequences from GenBank for alignments. Spaces or tabs are fine, but linebreaks are not, so remove those.",
    "hits graph": "The hits vs. e-value graph plots the number of hits from a BLAST search (total, number of genera, and number of species) at successively stricter e-value cutoffs, until hits from fewer than 10 species are found. This is useful for studying protein conservation and for determining what e-value to use for high quality alignments.",
    "align-o-tron": "BLASTer can automatically retrieve full-length protein sequences, then generate alignments and sequence logos for points on the hits vs. e-value graph with hits from fewer than 100 genera. Will only retrieve one hit per genus or, for points with fewer than 100 species, one hit per species.",
    "initial cutoff": "This is the e-value at which the hits vs. e-value graph and the align-o-tron will start counting down, expressed as 10^-X. Typically will be 0 or 1, but for genes with a very high number of homologs in GenBank, it may be helpful to start with a stricter e-value cutoff.",
    "number of hits": "BLAST searches have a maximum number of hits they can return. You would like this number to be higher than the number of total hits at the initial e-value cutoff. 20000 is a pretty safe bet, but particularly common proteins might return more hits than this.",
    "high quality alignment": "Takes a single e-value cutoff, but allows pruning results by species or genus, as well as various ways of trimming out sequences that are longer or shorter than desired.",
    "e-value cutoff": "Sets the stringency of the homology search. The hits vs. e-value graph is intended to give you the information you need to pick this.",
    "repeat step": "For several steps of the workflow, BLASTer gives the option to repeat the previous step with changed parameters, in case you’re not happy with the outcome.",
    "pruning": "BLASTer can pare down BLAST results, returning X hits per genus or species, as requested. The most common application is 1 hit per species, but 1 hit per genus can also be very useful.",
    "trimming": "BLASTer will report the length distribution of your BLAST results and give you the option to eliminate short or long sequences (or both), or to keep only long or short sequences.",
    "alignment": "MUSCLE is used to generate amino acid alignments, which are output in FASTA format and as sequence logos, for ease of visualization.",
    "homology table": "BLASTer can generate a table in which each homolog is BLASTed against the original sequence, reporting e-value, % identity, % similarity, and taxonomic information about the species from which that protein originates.",
    }

# FUNCTIONS

def helper(topic):
    """ Refers to help_dict for appropriate info, then allows user to check other topics, too. """
    print("\n-------------------------------------------\n")
    print("BLASTer Contextual Help System\n")
    answer = help_dict.get(topic, "That topic has not been implemented. Oops!")
    print(answer + "\n")

    topic_list = list(help_dict.keys())
    topics_string = ", ".join(topic_list)

    asking = True
    while asking:
        choice = input("Available help topics (or ‘X’ to exit):\n" + topics_string + "\n")
        if choice.upper() == "X" or choice.upper() == "EXIT" or choice.upper() == "":
            print("-------------------------------------------\n")
            asking = False
            return True

        answer = help_dict.get(choice, "Please select one of the listed topics.")
        print("\n" + answer + "\n")

def is_email(test):
    """ Checks a string to see if it looks like it might be an email address. """
    result = re.findall('\S+@\S+', test) # regular expression for any text flanking an "@"
    if len(result) == 0:
        return False
    else:
        return True

def remove_extension(filename):
    """Given a file name, returns the file name without an extension."""        
    parts = filename.split(".")
    assert len(parts) == 2
    return parts[0]

def get_species(title):
    """From a GenBank identifier line, extracts the full species and strain name."""
    in_name = 0
    species = ""
    if "MULTISPECIES" in title:
        species = "MULTISPECIES CONSENSUS"
    else:
        for character in title:
            if character in "[]":
                in_name = in_name + 1
            elif in_name == 1:
                species = species + character
            elif in_name >= 2:
                break
    return species

def species_name(identifier):  
    """From the output of get_species, extracts the binomial species name.""" 
    breakdown = identifier.split()     
    if len(breakdown) <= 1:    
        return identifier
    elif breakdown[0] == "Candidatus":          
        name = breakdown[1] + " " + breakdown[2]
        return name
    else:    
        name = breakdown[0] + " " + breakdown[1]      
        return name

def genus_name(identifier):   
    """From the output of get_species, extracts just the genus name.""" 
    breakdown = identifier.split()
    if len(breakdown) <= 1: 
        return identifier
    elif breakdown[0] == "Candidatus":      
        name = breakdown[1]
        return name
    else:
        name = breakdown[0]      
        return name

# CODE

# get the current date and reformat it
today = datetime.date.today()
dateyear = str(today.year)
if len(str(today.day)) == 2:
    dateday = str(today.day)
else:
    dateday = "0" + str(today.day)
datemo = str(today.month)
if len(str(today.month)) == 2:
    datemo = str(today.month)
else:
    datemo = "0" + str(today.month)
datestamp = dateyear + datemo + dateday

# print introductory header
print("BLASTer version %s" % version)
print(updated)
print("A tool for examining sequence conservation in proteins.")
print("by Michael J. Gray, M.S., Ph.D.\n")

# determine whether settings file exists
try:
    settings_handle = open("BLASTer_settings.txt", "r")
    
# if it doesn't, collect the necessary information
except: 
    print("No settings file has been created. You will need to enter")
    
    asking = True    
    while asking:
        email = input("1. a valid email address (required by NCBI): ")
        if email.upper() == "Q" or email.upper() == "QUIT":
            sys.exit()
        elif email.upper() == "H" or email.upper() == "HELP":
            helper("Entrez email")
            continue
        elif is_email(email) == False:
            print("This is not an email address.") # a quick, dumb test, but better than nothing
            sys.exit()
        else:
            asking == False
            break
              
    asking = True
    while asking:
        standalone_path = input("2. the path to the folder containing the standalone MUSCLE and blastp applications: ")
        if standalone_path.upper() == "Q" or standalone_path.upper() == "QUIT":
            sys.exit()
        elif standalone_path.upper() == "H" or standalone_path.upper() == "HELP":
            helper("standalones")
            continue
        elif os.path.exists(standalone_path) == False: # checks to see if this is an existing directory
            print("No such directory exists.")
            sys.exit()
        else:
            asking == False
            break
            
    asking = True
    while asking:
        basic_save_path = input("3. the path to the folder you want to save your data in: ")
        if standalone_path.upper() == "Q" or standalone_path.upper() == "QUIT":
            sys.exit()
        elif standalone_path.upper() == "H" or standalone_path.upper() == "HELP":
            helper("save path")
            continue
        elif os.path.exists(basic_save_path) == False: # checks to see if this is an existing directory
            print("No such directory exists.")
            sys.exit()
        else:
            asking == False
            break

    # write a new settings file
    outfile = open("BLASTer_settings.txt", "w")
    
    email_line = "email:" + email + "\n"
    standalone_line = "standalones:" + standalone_path + "\n"
    save_line = "basic_save_path:" + basic_save_path + "\n"
    
    outfile.write(email_line)
    outfile.write(standalone_line)
    outfile.write(save_line)

    outfile.close()

    print("\nFile created: BLASTer_settings.txt\n")

    settings_handle = open("BLASTer_settings.txt", "r")

# parse the settings file and tell the user what's in it

for line in settings_handle:
    line = line.rstrip() # get rid of the newline character
    separated = line.split(":") # break apart the line
    if separated[0] == "email":
        print("Entrez email = " + separated[1])
        Entrez.email = separated[1]
    elif separated[0] == "standalones":
        print("Standalone application directory = " + separated[1])
        standalone_path = separated[1]
    elif separated[0] == "basic_save_path":
        print("Output file directory = " + separated[1]) 
        basic_save_path = separated[1]
    else:
        print("Malformed settings file.") # this should never happen
        sys.exit()

# print helpful information and a separator to indicate that the actual program is starting

print("\nAt any prompt, enter 'Q' to quit or 'H' for contextual help.")
print("------------------------------------------------------------")

# clearly define where we are in the workflow and allow escape from nested loops
step = 1

while True:

# Step 1: Collect Initial Inputs
# inputs =
#   gene name = "name"
#   gene sequence = "input_record"

    while step == 1:

    # request the gene name
        asking = True
        while asking:
            name = input("Enter the protein name: ")
            if name.upper() == "Q" or name.upper() == "QUIT":
                sys.exit()
            elif name.upper() == "H" or name.upper() == "HELP":
                helper("protein name")
                continue
            else:
                asking = False

    # create a subdirectory to save the data in (and check to see if it already exists
        save_path = basic_save_path + "/" + datestamp + "_" + name + "/"
        try:
            os.makedirs(save_path)
        except OSError:
            asking = True
            while asking == True:
                question = input("A data file already exists. Do you want to overwrite it? Y or N: ")
                if question.upper() == "N" or question.upper() == "NO":
                    print("Please move or rename the existing file.")
                    sys.exit()
                elif question.upper() == "Y" or question.upper() == "YES":
                    asking = False
                    pass
                elif question.upper() == "Q" or question.upper() == "QUIT":
                    sys.exit()
                elif question.upper() == "H" or question.upper() == "HELP":
                    helper("overwriting data")
                    continue
                else:
                    continue

    # request the protein input sequence
        asking = True
        while asking:
            seq = input("Enter the full or partial protein sequence (with no linebreaks): ")
                # can't figure out how to get it to accept input with linebreaks, which worked with Python 2
            useq = seq.strip().replace(" ","").upper() # remove whitespace
            if useq == "Q" or useq == "QUIT":
                sys.exit()
            elif useq == "H" or useq == "HELP":
                helper("sequence")
                continue
            else:
                asking = False
                continue

    # convert the input sequence into a Biopython SeqRecord            
        input_seq = Seq(useq, IUPAC.protein)
        input_record = SeqRecord(input_seq, id=name, name=name, description="input sequence for BLAST search")

    # advance to the next step
        step = step + 1

# Step 2: The Align-O-Tron
#    inputs =
#       do you want a hits vs. evalue graph? (and associated tables)
#       if so:  initial e-value
#               upper limit of hits (recommend 20000)
#       do you want autoalignments and sequence logos?
#   outputs:
#       temporary XML file of raw BLAST results
#       requested tables and figures

    while step == 2:
        
    # do you want a hits vs. e-value graph?
        asking = True
        while asking: 
            question = input("\nGenerate a hits vs. e-value graph? Y or N: ")
            if question.upper() == "N" or question.upper() == "NO":
                alignotron = False
                asking = False
            elif question.upper() == "Y" or question.upper() == "YES":
                alignotron = True
                imagefilename = name + "_hits_graph.png"
                pathed_imagefilename = os.path.join(save_path, imagefilename)
                asking = False
            elif question.upper() == "Q" or question.upper() == "QUIT":
                sys.exit()
            elif question.upper() == "H" or question.upper() == "HELP":
                helper("hits graph")
                continue
            else:
                continue
            
    # do you want automatic alignments and sequence logos?
        if alignotron == True:
            asking = True
            while asking:
                autologos = input("Automatically generate alignments and logos? Y or N: ")
                if autologos.upper() == "N" or autologos.upper() == "NO":
                    want_alignments = False
                    asking = False
                elif autologos.upper() == "Y" or autologos.upper() == "YES":
                    want_alignments = True
                    asking = False
                elif autologos.upper() == "Q" or autologos.upper() == "QUIT":
                    sys.exit()
                elif autologos.upper() == "H" or autologos.upper() == "HELP":
                    helper("align-o-tron")
                    continue
                else:
                    continue

        # collect an initial e-value cutoff, and insist on an integer
            asking = True
            while asking: 
                try:
                    cutoff = input("Initial e-value cutoff (10^-X): ")
                    if str(cutoff).upper() == "Q" or str(cutoff).upper() == "QUIT":
                        sys.exit()
                    elif str(cutoff).upper() == "H" or str(cutoff).upper() == "HELP":
                        helper("initial cutoff")
                        continue
                    cutoff = int(cutoff) # in Python 3, input is always a string, and needs to be converted to a number
                except ValueError:
                    print("Please enter an integer.")
                    continue
                else:
                    asking = False

        # collect a value for upper limit on total hits, and insist on a positive integer
            asking = True
            while asking: 
                try:
                    limit = input("Maximum number of hits you would like returned (recommended = 20000): ")
                    if str(limit).upper() == "Q" or str(limit).upper() == "QUIT":
                        sys.exit()
                    elif str(limit).upper() == "H" or str(limit).upper() == "HELP":
                        helper("number of hits")
                        continue
                    elif str(limit).upper() == "":
                        limit = 20000
                    limit = int(limit) # in Python 3, input is always a string, and needs to be converted to a number
                except ValueError:
                    print("Please enter an integer.")
                    continue
                else:
                    if int(limit) < 0:
                        print("Please enter a positive number")
                        continue
                    else:
                        asking = False

        # initialize blank lists and tables for the align-o-tron

        # generate blank table of counts
            tablefile = name + "_alignotron_conservation_table.txt"
            pathed_tablefile = os.path.join(save_path, tablefile)
            save_file = open(pathed_tablefile, "w")
            first_row = name + " conservation table, generated " + datestamp + " using BLASTer version " + version + "\n"
            second_row = "Input sequence: " + useq + "\n"
            third_row = "E-value (10^-X) \t Total Hits \t Species \t Genera \n"
            save_file.write(first_row)
            save_file.write(second_row)
            save_file.write(third_row)
            save_file.close()

        # generate blank table of species
            species_tablefile = name + "_alignotron_species_table.txt"
            pathed_species_tablefile = os.path.join(save_path, species_tablefile)
            save_file = open(pathed_species_tablefile, "w")
            save_file.close()
            
        # generate blank table of genera
            genus_tablefile = name + "_alignotron_genus_table.txt"
            pathed_genus_tablefile = os.path.join(save_path, genus_tablefile)
            save_file = open(pathed_genus_tablefile, "w")
            save_file.close()

        # blank lists for figure generation
            evalue_axis = list()
            total_hits = list()
            species_hits = list()
            genus_hits = list()

        # set initial e-value cutoff (using input entered by user)
            e_thresh = math.pow(10, -int(cutoff))

        # generate initial BLAST result
            print("\nRunning BLAST search (this may take a while)...\n")
            search_handle = NCBIWWW.qblast("blastp", "nr", input_record.seq, expect=e_thresh, hitlist_size=limit)

        # write a temporary XML file containing the results
            xmlfilename = name + ".xml"
            pathed_xmlfilename = os.path.join(save_path, xmlfilename)
            xml_file = open(pathed_xmlfilename, "w")
            xml_file.write(search_handle.read())
            xml_file.close()
            search_handle.close()

        # begin homology detection loop
            tally = int((limit/2)*2) # cludgy way to set the value of "limit" to the value of "tally" while preserving "limit"
                # tally will update to reflect the number of species hits at each step of the loop    
            loop_counter = 0
            negative_log_e = cutoff - 1 # in each loop, will be +1, so this sets the first loop to the desired "cutoff"

            print("Beginning homology detection loop...\n")
            print(third_row.rstrip())
        
            while tally >= 10:       # stop script when results hit less than 10 species
                negative_log_e = negative_log_e + 1
                e_thresh = math.pow(10, -negative_log_e)
                loop_counter = loop_counter + 1

        # I have separated the data extraction into three steps below to simplify the logic of the code,
        # which otherwise would need a complicated and difficult to follow nested structure
        
            # count the number of total hits
                hit_count = 0
                search_handle = open(pathed_xmlfilename, "r")
                blast_record = SearchIO.read(search_handle, "blast-xml")
                for hit in blast_record:
                    for hsp in hit:
                        if hsp.evalue < e_thresh:
                            hit_count = hit_count + 1
                search_handle.close()
        
            # count the number of species hits and store their names       
                species_gi_list = list()
                species_counter = dict()
                search_handle = open(pathed_xmlfilename, "r")
                blast_record = SearchIO.read(search_handle, "blast-xml")
                for hit in blast_record:
                    for hsp in hit:
                        if hsp.evalue < e_thresh:
                            ident = get_species(hsp.hit_description)         
                            species = species_name(ident)
                            if ident.isspace() == False: # eliminates poorly annotated or misformatted hits
                                if ident != "MULTISPECIES CONSENSUS" and ident != "":     
                                    if species not in species_counter.keys():
                                        species_counter[species] = 1
                                        species_gi_list.append(hsp.hit_id)
                                    else:
                                        species_counter[species] = species_counter[species] + 1
                search_handle.close()
            
                tally = len(species_counter)       # see above - stop script when results hit less than 10 species

            # count the number of genus hits and store their names
                genus_gi_list = list()
                genus_counter = dict()
                search_handle = open(pathed_xmlfilename, "r")
                blast_record = SearchIO.read(search_handle, "blast-xml")
                for hit in blast_record:
                    for hsp in hit:
                        if hsp.evalue < e_thresh:
                            ident = get_species(hsp.hit_description)         
                            genus = genus_name(ident)
                            if ident.isspace() == False: # eliminates poorly annotated or misformatted hits
                                if ident != "MULTISPECIES CONSENSUS" and ident != "":     
                                    if genus not in genus_counter.keys():
                                        genus_counter[genus] = 1
                                        genus_gi_list.append(hsp.hit_id)
                                    else:
                                        genus_counter[genus] = genus_counter[genus] + 1
                search_handle.close()

            # generate data row for table of counts and print it
                row = str(negative_log_e) + "\t" + str(hit_count) + "\t" + str(len(species_counter)) + "\t" + str(len(genus_counter)) + "\n"
                save_file = open(pathed_tablefile, "a")    
                save_file.write(row)    
                save_file.close()
                print(row.rstrip())

            # add hit count data to graphing axes
                evalue_axis.append(negative_log_e)
                total_hits.append(hit_count)
                species_hits.append(len(species_counter))
                genus_hits.append(len(genus_counter))

            # generate table of species, overwriting previous version
                if loop_counter == 1:       # first loop - just write the initial file
                    save_file = open(pathed_species_tablefile, "w")
                    for species in species_counter:
                        row = species + "\t" + str(species_counter[species]) + "\n"
                        save_file.write(row)
                    save_file.close()
                
                else:       # subsequent loops - write secondary file, compare it to primary file, then write tertiary file
                    tempfile1 = os.path.join(save_path, "temp1.txt")
                    save_tempfile1 = open(tempfile1, "w")
                    for species in species_counter:
                        row = species + "\t" + str(species_counter[species]) + "\n"
                        save_tempfile1.write(row)
                    save_tempfile1.close()

                    old_list = open(pathed_species_tablefile, "r")
                    new_list = open(tempfile1, "r")
                    tempfile2 = os.path.join(save_path, "temp2.txt") 
                    save_file = open(tempfile2, "w")
                
                    changers = dict()
                    for line in new_list:
                        breakdown = line.rstrip("\n").split("\t")
                        changers[breakdown[0]] = breakdown[1]
                    new_list.close()        
            
                    for line in old_list:
                        words = line.rstrip("\n").split("\t")
                        first_word = words[0]
                        old_row = '\t'.join(words)
                        try:
                            new_number = changers[first_word]
                            row = old_row + "\t" + new_number + "\n"
                            save_file.write(row)
                        except KeyError:
                            row = old_row + "\t" + "0" + "\n"
                            save_file.write(row)  
                    save_file.close()
                    old_list.close()   
                             
                # overwrite old files 
                    os.remove(pathed_species_tablefile)
                    os.remove(tempfile1)
                    os.rename(tempfile2, pathed_species_tablefile)
    
            # generate table of genera, overwriting previous version
                if loop_counter == 1:       # first loop - just write the initial file
                    save_file = open(pathed_genus_tablefile, "w")
                    for genus in genus_counter:
                        row = genus + "\t" + str(genus_counter[genus]) + "\n"
                        save_file.write(row)
                    save_file.close()
                
                else:       # subsequent loops - write secondary file, compare it to primary file, then write tertiary file
                    tempfile1 = os.path.join(save_path, "temp1.txt")
                    save_tempfile1 = open(tempfile1, "w")
                    for genus in genus_counter:
                        row = genus + "\t" + str(genus_counter[genus]) + "\n"
                        save_tempfile1.write(row)
                    save_tempfile1.close()

                    old_list = open(pathed_genus_tablefile, "r")
                    new_list = open(tempfile1, "r")
                    tempfile2 = os.path.join(save_path, "temp2.txt") 
                    save_file = open(tempfile2, "w")
                
                    changers = dict()
                    for line in new_list:
                        breakdown = line.rstrip("\n").split("\t")
                        changers[breakdown[0]] = breakdown[1]
                    new_list.close()        
            
                    for line in old_list:
                        words = line.rstrip("\n").split("\t")
                        first_word = words[0]
                        old_row = '\t'.join(words)
                        try:
                            new_number = changers[first_word]
                            row = old_row + "\t" + new_number + "\n"
                            save_file.write(row)
                        except KeyError:
                            row = old_row + "\t" + "0" + "\n"
                            save_file.write(row)  
                    save_file.close()
                    old_list.close()   

                # overwrite old files 
                    os.remove(pathed_genus_tablefile)
                    os.remove(tempfile1)
                    os.rename(tempfile2, pathed_genus_tablefile)

            # decide whether to generate an alignment and request full-length sequences
                if want_alignments == True:

                # generating folders to store the data in                    
                    fasta_path = save_path + "/fasta_files/"    # make a folder for unaligned FASTA files
                    try:
                        os.makedirs(fasta_path)
                    except OSError:
                        pass

                    alignment_path = save_path + "/alignments/" # make a folder for FASTA alignments
                    try:
                        os.makedirs(alignment_path)
                    except OSError:
                        pass

                    logo_path = save_path + "/logos/"   # make a folder for sequence logos
                    try:
                        os.makedirs(logo_path)
                    except OSError:
                        pass

                # make a decision about whether to automatically generate an alignment       
                    decision = "no alignment"   # the default is NOT to make an alignment (because there are too many hits)

                    if len(genus_counter) < 100:    # align 1 per genus for large numbers of hits
                        decision = "genera"
                    if len(genus_counter) < 10 and len(species_counter) < 100:  # align 1 per species for small numbers of hits
                        decision = "species"

                    if decision != "no alignment":

                    # count hits and generate file names for FASTA and logo file, as well as titles for the logo images
                        if decision == "genera":
                            gi_str = ",".join(genus_gi_list)
                            fastafilename = name + "_e^-" + str(negative_log_e) + "_" + str(len(genus_counter)) + "genera.fasta"
                            graph_title = name + " alignment (" + str(len(genus_counter)) + " genera, e-value < 10^-" + str(negative_log_e) + ")"
                        elif decision == "species":
                            gi_str = ",".join(species_gi_list)
                            fastafilename = name + "_e^-" + str(negative_log_e) + "_" + str(len(species_counter)) + "species.fasta"
                            graph_title = name + " alignment (" + str(len(species_counter)) + " species, e-value < 10^-" + str(negative_log_e) + ")"
                            
                        else:
                            print("This should never happen!")
                            sys.exit()
                         
                        pathed_fastafilename = os.path.join(fasta_path, fastafilename)

                    # request full-length protein sequences from GenBank
                        handle = Entrez.efetch(db="protein", id=gi_str, rettype="fasta", retmode="text")        
                        records = SeqIO.parse(handle, "fasta")                                                 
                        SeqIO.write(records, pathed_fastafilename, "fasta") # writes a FASTA file with full-length protein sequences
                
                        if decision == "genera":
                            print(fastafilename + " saved, containing " + str(len(genus_gi_list)) + " " + decision + ".")
                        elif decision == "species":
                            print(fastafilename + " saved, containing " + str(len(species_gi_list)) + " " + decision + ".")
                        else:
                            print("This should never happen!")
                            sys.exit()

                    # generate actual alignment using MUSCLE
                        base_name = remove_extension(fastafilename)         
                        aln_name = base_name + ".aln"
                        pathed_aln_name = os.path.join(alignment_path, aln_name)
                        program = mstandalone_path + "/muscle3.8.31_i86darwin64"
                        
                        try:
                            cline = MuscleCommandline(program, input=pathed_fastafilename, out=pathed_aln_name)  
                            run = cline()
                        except:
                            print("Error! MUSCLE aligner not found in correct directory.")
                            sys.exit()
                        
                        print("MUSCLE alignment result file saved as %s" % aln_name)

                    # generate sequence logo using Logomaker and matplotlib
                        alignment_data = AlignIO.read(pathed_aln_name, "fasta") # parse the alignment file
                        
                        alignment_sequence_list = list()
                        for record in alignment_data:
                            alignment_sequence_list.append(str(record.seq)) # extract just the individual aligned sequences as strings                     

                        alignment_matrix = lm.alignment_to_matrix(sequences=alignment_sequence_list, to_type='counts', characters_to_ignore='.-X')

                    # short alignments look fine as is
                        if alignment_data.get_alignment_length() <= 60:
                            aln_logo = lm.Logo(alignment_matrix, color_scheme='chemistry', figsize=(12,2))
                            aln_logo.suptitle(graph_title, fontsize="14", fontweight="bold")
                            aln_logo.ax.set_ylabel('Counts')
                            aln_logo.style_spines(visible=False)
                            aln_logo.style_spines(spines=['left','bottom'], visible=True, linewidth=1)

                    # but logos with more than 60 stacks are unreadable, so need to be split into multiple rows
                        else: 

                        # make a temporary full-length logo to determine the appropriate y-axis limits
                            temp_logo = lm.Logo(alignment_matrix, color_scheme='chemistry', figsize=(12,2))
                            ylimits = temp_logo.ax.get_ylim()
                            plt.close()

                        # setting row height and alignment slicing start points
                            height_per_row = 2
                            
                            start_slice = 0
                            end_slice = 59
                            counter = 0
                            
                        # for alignments whose length is exactly divisible by 60
                            if alignment_data.get_alignment_length() % 60 == 0:
                                rows = alignment_data.get_alignment_length() // 60
                                total_height = height_per_row * rows
                                
                                fig = plt.figure(figsize=(12, total_height))
                                fig.suptitle(graph_title, fontsize="14", fontweight="bold")
                                spec = gridspec.GridSpec(ncols=1, nrows=rows, figure=fig)

                                for _ in range(rows):
                                    sub_ax = fig.add_subplot(spec[counter, :])
                                    aln_logo = lm.Logo(alignment_matrix.loc[start_slice:end_slice], ax=sub_ax, color_scheme='chemistry', figsize=(12,2))
                                    aln_logo.ax.set_ylim(ylimits)
                                    aln_logo.ax.set_ylabel('Counts')
                                    aln_logo.ax.set_anchor('W') # left-justify the plots
                                    aln_logo.style_spines(visible=False)
                                    aln_logo.style_spines(spines=['left','bottom'], visible=True, linewidth=1)
                                    aln_logo.style_xticks(anchor=start_slice, spacing=10)

                                    start_slice = end_slice + 1
                                    end_slice = end_slice + 60
                                    counter = counter + 1
                                    
                        # for alignments whose length is not exactly divisible by 60
                            else:
                                rows = (alignment_data.get_alignment_length() // 60) + 1   
                                ratio = (alignment_data.get_alignment_length() % 60) / 60
                                column1 = 12 * ratio
                                column2 = 12 - column1
                                total_height = height_per_row * rows

                                fig = plt.figure(figsize=(12, total_height))
                                fig.suptitle(graph_title, fontsize="14", fontweight="bold")
                                spec = gridspec.GridSpec(ncols=2, nrows=rows, width_ratios=[column1, column2], figure=fig)

                                for _ in range(rows-1):
                                    sub_ax = fig.add_subplot(spec[counter, :])
                                    aln_logo = lm.Logo(alignment_matrix.loc[start_slice:end_slice], ax=sub_ax, color_scheme='chemistry', figsize=(12,2))
                                    aln_logo.ax.set_ylim(ylimits)
                                    aln_logo.ax.set_ylabel('Counts')
                                    aln_logo.ax.set_anchor('W') # left-justify the plots
                                    aln_logo.style_spines(visible=False)
                                    aln_logo.style_spines(spines=['left','bottom'], visible=True, linewidth=1)
                                    aln_logo.style_xticks(anchor=start_slice, spacing=10)

                                    start_slice = end_slice + 1
                                    end_slice = end_slice + 60
                                    counter = counter + 1

                            # draw the short row
                                sub_ax = fig.add_subplot(spec[-1,0])
                                aln_logo = lm.Logo(alignment_matrix.loc[start_slice:], ax=sub_ax, color_scheme='chemistry', figsize=(column1,2))
                                aln_logo.ax.set_ylim(ylimits)
                                aln_logo.ax.set_ylabel('Counts')
                                aln_logo.ax.set_anchor('W') # left-justify the plots
                                aln_logo.style_spines(visible=False)
                                aln_logo.style_spines(spines=['left','bottom'], visible=True, linewidth=1)
                                aln_logo.style_xticks(anchor=start_slice, spacing=10)

                    # save logo
                        logo_name = base_name + "_logo.png"
                        pathed_logo_name = os.path.join(logo_path, logo_name)
                        plt.savefig(pathed_logo_name, dpi=300)
                        plt.close()

                        print("Sequence logo saved as %s" % logo_name)

        # draw the hits vs. e-value graph with matplotlib
            plt.figure(1)
            plt.subplot(111)
            plt.title('BLAST hits vs. E-Value')
            plt.semilogy(evalue_axis, total_hits, 'k-', linewidth=2.0, label='Total Hits')
            plt.semilogy(evalue_axis, species_hits, 'b-', linewidth=2.0, label='Species')
            plt.semilogy(evalue_axis, genus_hits, 'r-', linewidth=2.0, label='Genera')
            plt.legend(loc='upper right')
            plt.xlabel(r'E-Value (10$^\mathrm{-X}$)')  # renders the "-X" as a superscript
            plt.ylabel('BLAST hits')
            
        # save graph as a PNG image
            plt.savefig(pathed_imagefilename, dpi=300)
            plt.close()
            print ("\nHits vs. E-Value graph saved as " + imagefilename)
            
        # get rid of XML file
            os.remove(pathed_xmlfilename) 
        
    # do you want a high quality alignment?
        asking = True
        while asking: 
            question = input("\nGenerate a high quality alignment? Y or N: ")
            if question.upper() == "N" or question.upper() == "NO":
                sys.exit()
            elif question.upper() == "Y" or question.upper() == "YES":
                step = step + 1
                asking = False
            elif question.upper() == "Q" or question.upper() == "QUIT":
                sys.exit()
            elif question.upper() == "H" or question.upper() == "HELP":
                helper("high quality alignment")
                continue
            else:
                continue

# Step 3: 
#    inputs =
#       e-value threshold
#       upper limit of hits 
#   outputs:
#       temporary XML file of raw BLAST results
#       number of total hits
#       tables of species and genus frequency

    while step == 3:

    # collect an e-value cutoff, and insist on an integer
        asking = True
        while asking: 
            try:
                cutoff = input("Enter an e-value cutoff (10^-X): ")
                if str(cutoff).upper() == "Q" or str(cutoff).upper() == "QUIT":
                    sys.exit()
                elif str(cutoff).upper() == "H" or str(cutoff).upper() == "HELP":
                    helper("e-value cutoff")
                    continue
                cutoff = int(cutoff) # in Python 3, input is always a string, and needs to be converted to a number
            except ValueError:
                print("Please enter an integer.")
                continue
            else:
                asking = False

    # collect an number of hits, and insist on an integer
        asking = True
        while asking: 
            try:
                hits = input("How many hits would you like you returned? ")
                if str(hits).upper() == "Q" or str(hits).upper() == "QUIT":
                    sys.exit()
                elif str(hits).upper() == "H" or str(hits).upper() == "HELP":
                    helper("number of hits")
                    continue
                hits = int(hits) # in Python 3, input is always a string, and needs to be converted to a number
            except ValueError:
                print("Please enter an integer.")
                continue
            else:
                asking = False

    # set initial e-value cutoff (using input entered by user)
        e_thresh = math.pow(10, -int(cutoff))

    # generate initial BLAST result
        print("\nRunning BLAST search (this may take a while)...\n")
        search_handle = NCBIWWW.qblast("blastp", "nr", input_record.seq, expect=e_thresh, hitlist_size=hits)

    # write a temporary XML file containing the results
        xmlfilename = name + ".xml"
        pathed_xmlfilename = os.path.join(save_path, xmlfilename)
        xml_file = open(pathed_xmlfilename, "w")
        xml_file.write(search_handle.read())
        xml_file.close()
        search_handle.close()

    # count the number of total hits
        hit_count = 0
        search_handle = open(pathed_xmlfilename, "r")
        blast_record = SearchIO.read(search_handle, "blast-xml")
        for hit in blast_record:
            hit_count = hit_count + 1
        print("Results: " + str(hit_count) + " hits.")

    # make table of species frequency
        tablefile = name + "_species_table_" + datestamp + ".txt" 
        pathed_tablefile = os.path.join(save_path, tablefile) 
        save_file = open(pathed_tablefile, "w") 
        first_row = name + " species frequency table, generated " + datestamp + " using BLASTer version " + version + "\n"
        second_row = "Input sequence: " + useq + "\n"  
        third_row = "E-Value less than: " + str(e_thresh) + "\n"     
        fourth_row = "Total Hits: " + str(hit_count) + "\n"     
        fifth_row = "Species \t Occurance\n"
        save_file.write(first_row)    
        save_file.write(second_row)   
        save_file.write(third_row)   
        save_file.write(fourth_row)   
        save_file.write(fifth_row)   

        species_counter = dict()
        search_handle = open(pathed_xmlfilename, "r")
        blast_record = SearchIO.read(search_handle, "blast-xml")
        for hit in blast_record:
            for hsp in hit:
                ident = get_species(hsp.hit_description)         
                species = species_name(ident)
                if ident.isspace() == False: # eliminates poorly annotated or misformatted hits
                    if ident != "MULTISPECIES CONSENSUS" and ident != "":     
                        if species not in species_counter.keys():
                            species_counter[species] = 1
                        else:
                            species_counter[species] = species_counter[species] + 1

        for species in species_counter:
            row = species + "\t" + str(species_counter[species]) + "\n"
            save_file.write(row)

        print("Species occurance table for " + name + " at e-value < " + str(e_thresh) + " (" + str(len(species_counter)) + " species) saved as " + tablefile)
        save_file.close()
        search_handle.close()
        
    # make table of genus frequency
        tablefile = name + "_genus_table_" + datestamp + ".txt" 
        pathed_tablefile = os.path.join(save_path, tablefile) 
        save_file = open(pathed_tablefile, "w") 
        first_row = name + " genus frequency table, generated " + datestamp + " using BLASTer version " + version + "\n"
        second_row = "Input sequence: " + useq + "\n"  
        third_row = "E-Value less than: " + str(e_thresh) + "\n"     
        fourth_row = "Total Hits: " + str(hit_count) + "\n"     
        fifth_row = "Genus \t Occurance \n"
        save_file.write(first_row)    
        save_file.write(second_row)   
        save_file.write(third_row)   
        save_file.write(fourth_row)   
        save_file.write(fifth_row)   

        genus_counter = dict()
        search_handle = open(pathed_xmlfilename, "r")
        blast_record = SearchIO.read(search_handle, "blast-xml")
        for hit in blast_record:
            for hsp in hit:
                ident = get_species(hsp.hit_description)         
                genus = genus_name(ident)
                if ident.isspace() == False: # eliminates poorly annotated or misformatted hits
                    if ident != "MULTISPECIES CONSENSUS" and ident != "":     
                        if genus not in genus_counter.keys():
                            genus_counter[genus] = 1
                        else:
                            genus_counter[genus] = genus_counter[genus] + 1
                            
        for genus in genus_counter:
            row = genus + "\t" + str(genus_counter[genus]) + "\n"
            save_file.write(row)

        print("Genus occurance table for " + name + " at e-value < " + str(e_thresh) + " (" + str(len(genus_counter)) + " genera) saved as " + tablefile)
        save_file.close()
        search_handle.close()

    # request clarification about what to do next
        asking = True
        while asking:
            clarify = input("\n(C)ontinue to data extraction step or (R)epeat search with changed parameters? ") 
            if clarify.upper() == "C" or clarify.upper() == "CONTINUE":
                step = step + 1
                break
            elif clarify.upper() == "R" or clarify.upper() == "REPEAT":
                asking = False 
            elif clarify.upper() == "H" or clarify.upper() == "HELP":
                helper("repeat step")
                continue
            elif clarify.upper() == "Q" or clarify.upper() == "QUIT":
                sys.exit()
            else:
                continue

# Step 4: Pruning the List & Converting the Output File
# inputs =
#   pruning?
#   species or genus?
#   how many per species / genus?
# outputs =
#   FASTA formatted full-length protein list file
#   number of proteins

    while step == 4:
        gi_list = list()

        asking = True
        while asking:
            query = input("\nWould you like to generate a pruned result file? (Y/N): ")
            if query.upper() == "Q" or query.upper() == "QUIT":
                sys.exit()
            if query.upper() == "H" or query.upper() == "HELP":
                helper("pruning")
                continue

         # if no pruning is requested, just return every hit in the XML file
            elif query.upper() == "N":  
                print("Extracting relevant data...")
                search_handle = open(pathed_xmlfilename, "r")
                fastafilename = name + "_" + datestamp + ".fasta" 
                pathed_fastafilename = os.path.join(save_path, fastafilename) 
                blast_record = SearchIO.read(search_handle, "blast-xml")
                for hit in blast_record:
                    for hsp in hit:
                        gi_list.append(hsp.hit_id)
                search_handle.close()
                asking = False

        # if pruning is requested, determine how that pruning should be done
            elif query.upper() == "Y":
                species_counter = dict()

                asking2 = True
                while asking2:
                    query2 = input("Prune by species or by genus? (S/G): ")

                    if query2.upper() == "Q" or query2.upper() == "QUIT":
                        sys.exit()
                    elif query2.upper() == "H" or query2.upper() == "HELP":
                        helper("pruning")
                        continue
            
                    elif query2.upper() == "S" or query2.upper() == "SPECIES": # extract the first "number" hits per species from the XML file
                        asking3 = True
                        while asking3:
                            try:
                                number = input("How many hits per species? ")
                                if str(number).upper() == "Q" or str(number).upper() == "QUIT":
                                    sys.exit()
                                elif str(number).upper() == "H" or str(number).upper() == "HELP":
                                    helper("pruning")
                                    continue
                                number = int(number)
                            except ValueError:
                                print("Please enter an integer.")
                                continue
                            else:
                                asking3 = False
                
                        print("Extracting relevant data...")
                        search_handle = open(pathed_xmlfilename, "r")
                        fastafilename = name + "_" + datestamp + ".fasta" 
                        pathed_fastafilename = os.path.join(save_path, fastafilename) 
                        blast_record = SearchIO.read(search_handle, "blast-xml")
                        for hit in blast_record:
                            for hsp in hit:
                                ident = get_species(hsp.hit_description)         
                                species = species_name(ident)
                                if ident.isspace() == False: # eliminates poorly annotated or misformatted hits
                                    if ident != "MULTISPECIES CONSENSUS" and ident != "":     
                                        if species not in species_counter.keys():
                                            species_counter[species] = 1
                                            gi_list.append(hsp.hit_id)
                                        else:
                                            if species_counter[species] < number:
                                                gi_list.append(hsp.hit_id)
                                                species_counter[species] = species_counter[species] + 1
                                            else:
                                                species_counter[species] = species_counter[species] + 1
                        search_handle.close()
                        asking2 = False
                        asking = False

                    elif query2.upper() == "G" or query2.upper() == "GENUS": # extract the first "number" hits per genus from the XML file
                        asking3 = True
                        while asking3:
                            try:
                                number = input("How many hits per genus? ")
                                if str(number).upper() == "Q" or str(number).upper() == "QUIT":
                                    sys.exit()
                                elif str(number).upper() == "H" or str(number).upper() == "HELP":
                                    helper("pruning")
                                    continue
                                number = int(number)
                            except ValueError:
                                print("Please enter an integer.")
                                continue
                            else:
                                asking3 = False
                
                        print("Extracting relevant data...")
                        search_handle = open(pathed_xmlfilename, "r")
                        fastafilename = name + "_" + datestamp + ".fasta" 
                        pathed_fastafilename = os.path.join(save_path, fastafilename) 
                        blast_record = SearchIO.read(search_handle, "blast-xml")
                        for hit in blast_record:
                            for hsp in hit:
                                ident = get_species(hsp.hit_description)         
                                species = genus_name(ident)
                                if ident.isspace() == False: # eliminates poorly annotated or misformatted hits
                                    if ident != "MULTISPECIES CONSENSUS" and ident != "":     
                                        if species not in species_counter.keys():
                                            species_counter[species] = 1
                                            gi_list.append(hsp.hit_id)
                                        else:
                                            if species_counter[species] < number:
                                                gi_list.append(hsp.hit_id)
                                                species_counter[species] = species_counter[species] + 1
                                            else:
                                                species_counter[species] = species_counter[species] + 1
                        search_handle.close()
                        asking2 = False
                        asking = False

                    else:
                        continue

    # request full-length protein sequences from GenBank
        print("Requesting GenBank files...")                                                 
        gi_str = ",".join(gi_list)                                                          
        handle = Entrez.efetch(db="protein", id=gi_str, rettype="fasta", retmode="text")        

    # and save them into a FASTA file
        print("Reformatting results to FASTA file, containing " + str(len(gi_list)) + " records: %s" % fastafilename)
        records = SeqIO.parse(handle, "fasta")                                                 
        SeqIO.write(records, pathed_fastafilename, "fasta")
            
    # request clarification about what to do next
        asking = True
        while asking:
            clarify = input("\n(C)ontinue to trimming step or (R)epeat pruning with changed parameters? ")
            if clarify.upper() == "C" or clarify.upper() == "CONTINUE":
                step = step + 1
                os.remove(pathed_xmlfilename)
                break
            elif clarify.upper() == "R" or clarify.upper() == "REPEAT":
                asking = False 
            elif clarify.upper() == "Q" or clarify.upper() == "QUIT":
                os.remove(pathed_xmlfilename)
                sys.exit()
            elif clarify.upper() == "H" or clarify.upper() == "HELP":
                helper("repeat step")
                continue
            else:
                continue
            
# Step 5: Trim Contents of Protein List By Length?
#
# outputs = 
#   length distribution (relative to full length of top hit)
# inputs = 
#   trim or continue to next step?
#   cutoff length method
# outputs = 
#   trimmed FASTA formatted protein list file (overwrites untrimmed file)
#   number of protins

    while step == 5:

    # determine the length of the top hit (which will probably be identical to the input sequence)
        first_record = next(SeqIO.parse(pathed_fastafilename, "fasta"))
        length = len(first_record)
        print("\n" + name + " is a " + str(length) + " amino acid protein. " + fastafilename + " contains " + str(len(gi_list)) + " sequences,")
       
    # determine the length distribution of the sample
        distribution = list()
        for seq_record in SeqIO.parse(pathed_fastafilename, "fasta"):
            seq = seq_record.seq
            distribution.append(len(seq))
        print("with an average length of " + str(np.mean(distribution)) + " and standard deviation of " + str(np.std(distribution)) + ".")

    # calculate basic percentage cutoff lengths
        short_cutoff = (length * 75) / 100
        long_cutoff = (length * 125) / 100
        
    # count sequences and report % long and short sequences
        short_count = 0
        long_count = 0
        for seq_record in SeqIO.parse(pathed_fastafilename, "fasta"):
            seq = seq_record.seq
            if len(seq) < short_cutoff:
                short_count = short_count + 1
            elif len(seq) > long_cutoff:
                long_count = long_count + 1                
        print(str(short_count) + " are sequences less than " + str(short_cutoff) + " amino acids (75% of " + name + ")")
        print(str(long_count) + " are sequences more than " + str(long_cutoff) + " amino acids (125% of " + name + ")\n")

    # request clarification about what to do next, and trim if required
        asking = True
        while asking:
            clarify = input("Choose one:\nTrim (S)hort sequences\nTrim (L)ong sequences\nTrim (B)oth long and short sequences\nKeep only long sequences (OL)\nKeep only short sequences (OS)\nDo (N)o trimming\n")
            
            if clarify.upper() == "Q" or clarify.upper() == "QUIT":
                sys.exit()
            elif clarify.upper() == "H" or clarify.upper() == "HELP":
                helper("trimming")
                continue

        # choice == eliminate short sequences
            elif clarify.upper() == "S" or clarify.upper() == "SHORT":  
                
            # ask about how to trim (numerically)
                asking2 = True
                while asking2:
                    ask = input("Trim by (P)ercent length (75%), (S)tandard deviation, or (D)efined length? ")
                    if ask.upper() == "Q" or ask.upper() == "QUIT":
                        sys.exit()
                    elif ask.upper() == "H" or ask.upper() == "HELP":
                        helper("trimming")
                        continue
                    elif ask.upper() == "P" or ask.upper() == "PERCENT" or ask.upper() == "PERCENT LENGTH":
                        asking2 = False # leave short_cutoff as 75%
                    elif ask.upper() == "S" or ask.upper() == "SD" or ask.upper() == "STANDARD DEVIATION":
                        asking3 = True
                        while asking3:
                            try:
                                SD = input("How many standard deviations? ")
                                if str(SD).upper() == "Q" or str(SD).upper() == "QUIT":
                                    sys.exit()
                                elif str(SD).upper() == "H" or str(SD).upper() == "HELP":
                                    helper("trimming")
                                    continue
                                SD = int(SD)
                            except ValueError:
                                print("Please enter an integer.")
                                continue
                            else:
                                asking3 = False
                        short_cutoff = np.mean(distribution) - (int(SD) * np.std(distribution))
                        asking2 = False
                    elif ask.upper() == "D" or ask.upper() == "DEFINED" or ask.upper() == "DEFINED LENGTH":
                        asking3 = True
                        while asking3:
                            try:
                                cutoff1 = input("Minimum length? ")
                                if str(cutoff1).upper() == "Q" or str(cutoff1).upper() == "QUIT":
                                    sys.exit()
                                elif str(cutoff1).upper() == "H" or str(cutoff1).upper() == "HELP":
                                    helper("trimming")
                                    continue
                                cutoff1 = int(cutoff1)
                            except ValueError:
                                print("Please enter an integer.")
                                continue
                            else:
                                asking3 = False
                        short_cutoff = int(cutoff1) - 1 # make the output more intuitive
                        asking2 = False
                    else:
                        continue
                    
            # trim short sequences
                search_handle = open(pathed_fastafilename, "r")
                pathed_tempfilename = os.path.join(save_path, "temp.fasta") 
                savefile = open(pathed_tempfilename, "w")
                counter = 0
                for seq_record in SeqIO.parse(pathed_fastafilename, "fasta"):
                    seq = seq_record.seq
                    if len(seq) > short_cutoff:
                        counter = counter + 1
                        SeqIO.write(seq_record, savefile, "fasta")
                savefile.close()
                search_handle.close()
            # overwrite old file 
                os.remove(pathed_fastafilename)
                os.rename(pathed_tempfilename, pathed_fastafilename)
            # print file name and new length
                print(fastafilename + " trimmed to " + str(counter) + " sequences longer than " + str(short_cutoff) + " amino acids.")
                step = step + 1
                break
                
        # choice == eliminate long sequences
            elif clarify.upper() == "L" or clarify.upper() == "LONG":
                
            # ask about how to trim
                asking2 = True
                while asking2:
                    ask = input("Trim by (P)ercent length (125%), (S)tandard deviation, or (D)efined length? ")
                    if ask.upper() == "Q" or ask.upper() == "QUIT":
                        sys.exit()
                    elif ask.upper() == "H" or ask.upper() == "HELP":
                        helper("trimming")
                        continue
                    elif ask.upper() == "P" or ask.upper() == "PERCENT" or ask.upper() == "PERCENT LENGTH":
                        asking2 = False # leave long_cutoff as 125%
                    elif ask.upper() == "S" or ask.upper() == "SD" or ask.upper() == "STANDARD DEVIATION":
                        asking3 = True
                        while asking3:
                            try:
                                SD = input("How many standard deviations? ")
                                if str(SD).upper() == "Q" or str(SD).upper() == "QUIT":
                                    sys.exit()
                                elif str(SD).upper() == "H" or str(SD).upper() == "HELP":
                                    helper("trimming")
                                    continue
                                SD = int(SD)
                            except ValueError:
                                print("Please enter an integer.")
                                continue
                            else:
                                asking3 = False
                        long_cutoff = np.mean(distribution) + (int(SD) * np.std(distribution))
                        asking2 = False
                    elif ask.upper() == "D" or ask.upper() == "DEFINED" or ask.upper() == "DEFINED LENGTH":
                        asking3 = True
                        while asking3:
                            try:
                                cutoff1 = input("Maximum length? ")
                                if str(cutoff1).upper() == "Q" or str(cutoff1).upper() == "QUIT":
                                    sys.exit()
                                elif str(cutoff1).upper() == "H" or str(cutoff1).upper() == "HELP":
                                    helper("trimming")
                                    continue
                                cutoff1 = int(cutoff1)
                            except ValueError:
                                print("Please enter an integer.")
                                continue
                            else:
                                asking3 = False
                        short_cutoff = int(cutoff1) + 1 # make the output more intuitive
                        asking2 = False
                    else:
                        continue
                        
             # trim long sequences
                search_handle = open(pathed_fastafilename, "r")
                pathed_tempfilename = os.path.join(save_path, "temp.fasta") 
                savefile = open(pathed_tempfilename, "w")
                counter = 0
                for seq_record in SeqIO.parse(pathed_fastafilename, "fasta"):
                    seq = seq_record.seq
                    if len(seq) < long_cutoff:
                        counter = counter + 1
                        SeqIO.write(seq_record, savefile, "fasta")
                savefile.close()
                search_handle.close()
                
            # overwrite old file 
                os.remove(pathed_fastafilename)
                os.rename(pathed_tempfilename, pathed_fastafilename)
                
            # print file name and new length
                print(fastafilename + " trimmed to " + str(counter) + " sequences shorter than " + str(long_cutoff) + " amino acids.")
                step = step + 1
                break
            
        # choice == eliminate both short and long sequences        
            elif clarify.upper() == "B" or clarify.upper() == "BOTH":
                
            # ask about how to trim
                asking2 = True

                asking2 = True
                while asking2:
                    ask = input("Trim by (P)ercent length (75-125%), (S)tandard deviation, or (D)efined length? ")
                    if ask.upper() == "Q" or ask.upper() == "QUIT":
                        sys.exit()
                    elif ask.upper() == "H" or ask.upper() == "HELP":
                        helper("trimming")
                        continue
                    elif ask.upper() == "P" or ask.upper() == "PERCENT" or ask.upper() == "PERCENT LENGTH":
                        asking2 = False # leave short_cutoff as 75% and long_cutoff as 125%
                    elif ask.upper() == "S" or ask.upper() == "SD" or ask.upper() == "STANDARD DEVIATION":
                        asking3 = True
                        while asking3:
                            try:
                                SD = input("How many standard deviations? ")
                                if str(SD).upper() == "Q" or str(SD).upper() == "QUIT":
                                    sys.exit()
                                elif str(SD).upper() == "H" or str(SD).upper() == "HELP":
                                    helper("trimming")
                                    continue
                                SD = int(SD)
                            except ValueError:
                                print("Please enter an integer.")
                                continue
                            else:
                                asking3 = False
                        short_cutoff = np.mean(distribution) - (int(SD) * np.std(distribution))
                        long_cutoff = np.mean(distribution) + (int(SD) * np.std(distribution))
                        asking2 = False
                    elif ask.upper() == "D" or ask.upper() == "DEFINED" or ask.upper() == "DEFINED LENGTH":
                        asking3 = True
                        while asking3:
                            try:
                                cutoff1 = input("Minimum length? ")
                                if str(cutoff1).upper() == "Q" or str(cutoff1).upper() == "QUIT":
                                    sys.exit()
                                elif str(cutoff1).upper() == "H" or str(cutoff1).upper() == "HELP":
                                    helper("trimming")
                                    continue
                                cutoff1 = int(cutoff1)
                            except ValueError:
                                print("Please enter an integer.")
                                continue
                            else:
                                asking3 = False
                        short_cutoff = int(cutoff1) - 1 # make the output more intuitive
                        asking4 = True
                        while asking4:
                            try:
                                cutoff2 = input("Maximum length? ")
                                if str(cutoff2).upper() == "Q" or str(cutoff2).upper() == "QUIT":
                                    sys.exit()
                                elif str(cutoff2).upper() == "H" or str(cutoff2).upper() == "HELP":
                                    helper("trimming")
                                    continue
                                cutoff2 = int(cutoff2)
                            except ValueError:
                                print("Please enter an integer.")
                                continue
                            else:
                                asking4 = False
                        short_cutoff = int(cutoff1) + 1 # make the output more intuitive
                        asking2 = False
                    else:
                        continue

            # trim both long and short sequences
                search_handle = open(pathed_fastafilename, "r")
                pathed_tempfilename = os.path.join(save_path, "temp.fasta") 
                savefile = open(pathed_tempfilename, "w")
                counter = 0
                for seq_record in SeqIO.parse(pathed_fastafilename, "fasta"):
                    seq = seq_record.seq
                    if len(seq) > short_cutoff and len(seq) < long_cutoff:
                        counter = counter + 1
                        SeqIO.write(seq_record, savefile, "fasta")
                savefile.close()
                search_handle.close()
                
            # overwrite old file 
                os.remove(pathed_fastafilename)
                os.rename(pathed_tempfilename, pathed_fastafilename)
                
            # print file name and new length
                print(fastafilename + " trimmed to " + str(counter) + " sequences between " + str(short_cutoff) + " and " + str(long_cutoff) + " amino acids.")
                step = step + 1
                break

        # choice == keep only long sequences
            elif clarify.upper() == "OL" or clarify.upper() == "ONLY LONG":
                
            # ask about how to trim (numerically)
                asking2 = True
                while asking2:
                    ask = input("Trim by (P)ercent length (>125%), (S)tandard deviation, or (D)efined length? ")
                    if ask.upper() == "Q" or ask.upper() == "QUIT":
                        sys.exit()
                    elif ask.upper() == "H" or ask.upper() == "HELP":
                        helper("trimming")
                        continue
                    elif ask.upper() == "P" or ask.upper() == "PERCENT" or ask.upper() == "PERCENT LENGTH":
                        short_cutoff = (length * 125) / 100
                        asking2 = False 
                    elif ask.upper() == "S" or ask.upper() == "SD" or ask.upper() == "STANDARD DEVIATION":
                        asking3 = True
                        while asking3:
                            try:
                                SD = input("How many standard deviations? ")
                                if str(SD).upper() == "Q" or str(SD).upper() == "QUIT":
                                    sys.exit()
                                elif str(SD).upper() == "H" or str(SD).upper() == "HELP":
                                    helper("trimming")
                                    continue
                                SD = int(SD)
                            except ValueError:
                                print("Please enter an integer.")
                                continue
                            else:
                                asking3 = False
                        short_cutoff = np.mean(distribution) + (int(SD) * np.std(distribution))
                        asking2 = False
                    elif ask.upper() == "D" or ask.upper() == "DEFINED" or ask.upper() == "DEFINED LENGTH":
                        asking3 = True
                        while asking3:
                            try:
                                cutoff1 = input("Minimum length? ")
                                if str(cutoff1).upper() == "Q" or str(cutoff1).upper() == "QUIT":
                                    sys.exit()
                                elif str(cutoff1).upper() == "H" or str(cutoff1).upper() == "HELP":
                                    helper("trimming")
                                    continue
                                cutoff1 = int(cutoff1)
                            except ValueError:
                                print("Please enter an integer.")
                                continue
                            else:
                                asking3 = False
                        short_cutoff = int(cutoff1) + 1 # make the output more intuitive
                        asking2 = False
                    else:
                        continue
                    
            # trim short sequences
                search_handle = open(pathed_fastafilename, "r")
                pathed_tempfilename = os.path.join(save_path, "temp.fasta") 
                savefile = open(pathed_tempfilename, "w")
                counter = 0
                for seq_record in SeqIO.parse(pathed_fastafilename, "fasta"):
                    seq = seq_record.seq
                    if len(seq) > short_cutoff:
                        counter = counter + 1
                        SeqIO.write(seq_record, savefile, "fasta")
                savefile.close()
                search_handle.close()
            # overwrite old file 
                os.remove(pathed_fastafilename)
                os.rename(pathed_tempfilename, pathed_fastafilename)
            # print file name and new length
                print(fastafilename + " trimmed to " + str(counter) + " sequences longer than " + str(short_cutoff) + " amino acids.")
                step = step + 1
                break
                 
        # choice == keep only short sequences
            elif clarify.upper() == "OS" or clarify.upper() == "ONLY SHORT":
                
            # ask about how to trim
                asking2 = True
                while asking2:
                    ask = input("Trim by (P)ercent length (<75%), (S)tandard deviation, or (D)efined length? ")
                    if ask.upper() == "Q" or ask.upper() == "QUIT":
                        sys.exit()
                    elif ask.upper() == "H" or ask.upper() == "HELP":
                        helper("trimming")
                        continue
                    elif ask.upper() == "P" or ask.upper() == "PERCENT" or ask.upper() == "PERCENT LENGTH":
                        long_cutoff = (length * 75) / 100
                        asking2 = False 
                    elif ask.upper() == "S" or ask.upper() == "SD" or ask.upper() == "STANDARD DEVIATION":
                        asking3 = True
                        while asking3:
                            try:
                                SD = input("How many standard deviations? ")
                                if str(SD).upper() == "Q" or str(SD).upper() == "QUIT":
                                    sys.exit()
                                elif str(SD).upper() == "H" or str(SD).upper() == "HELP":
                                    helper("trimming")
                                    continue
                                SD = int(SD)
                            except ValueError:
                                print("Please enter an integer.")
                                continue
                            else:
                                asking3 = False
                        long_cutoff = np.mean(distribution) - (int(SD) * np.std(distribution))
                        asking2 = False
                    elif ask.upper() == "D" or ask.upper() == "DEFINED" or ask.upper() == "DEFINED LENGTH":
                        asking3 = True
                        while asking3:
                            try:
                                cutoff1 = input("Maximum length? ")
                                if str(cutoff1).upper() == "Q" or str(cutoff1).upper() == "QUIT":
                                    sys.exit()
                                elif str(cutoff1).upper() == "H" or str(cutoff1).upper() == "HELP":
                                    helper("trimming")
                                    continue
                                cutoff1 = int(cutoff1)
                            except ValueError:
                                print("Please enter an integer.")
                                continue
                            else:
                                asking3 = False
                        short_cutoff = int(cutoff1) - 1 # make the output more intuitive
                        asking2 = False
                    else:
                        continue
                        
             # trim long sequences
                search_handle = open(pathed_fastafilename, "r")
                pathed_tempfilename = os.path.join(save_path, "temp.fasta") 
                savefile = open(pathed_tempfilename, "w")
                counter = 0
                for seq_record in SeqIO.parse(pathed_fastafilename, "fasta"):
                    seq = seq_record.seq
                    if len(seq) < long_cutoff:
                        counter = counter + 1
                        SeqIO.write(seq_record, savefile, "fasta")
                savefile.close()
                search_handle.close()
                
            # overwrite old file 
                os.remove(pathed_fastafilename)
                os.rename(pathed_tempfilename, pathed_fastafilename)
                
            # print file name and new length
                print(fastafilename + " trimmed to " + str(counter) + " sequences shorter than " + str(long_cutoff) + " amino acids.")
                step = step + 1
                break

            elif clarify.upper() == "N" or clarify.upper() == "NO" or clarify.upper() == "NONE":
                step = step + 1
                break

            else:
                continue

# Step 6: Generate Alignment
# inputs = 
#   quit, skip step, or continue?
#   FASTA formatted gene list file from previous step
# outputs = 
#   MUSCLE-aligned, FASTA-formatted alignment file
#   sequence logo

    while step == 6:
        
    # request clarification about what to do next
        asking = True
        while asking:
            clarify = input("\n(C)ontinue to alignment step or (S)kip alignment? ")
            
            if clarify.upper() == "Q" or clarify.upper() == "QUIT":
                sys.exit()
            elif clarify.upper() == "H" or clarify.upper() == "HELP":
                helper("alignment")
                continue
            elif clarify.upper() == "S" or clarify.upper() == "SKIP":
                step = step + 1
                break
            
            elif clarify.upper() == "C" or clarify.upper() == "CONTINUE":

            # generate actual alignment using MUSCLE
                base_name = remove_extension(fastafilename)         
                aln_name = base_name + ".aln"
                pathed_aln_name = os.path.join(save_path, aln_name)
                program = standalone_path + "/muscle3.8.31_i86darwin64"
                        
                try:
                    cline = MuscleCommandline(program, input=pathed_fastafilename, out=pathed_aln_name)  
                    run = cline()
                except:
                    print("Error! MUSCLE aligner not found in correct directory.")
                    sys.exit()
                        
                print("MUSCLE alignment result file saved as %s" % aln_name)

            # generate sequence logo using Logomaker and matplotlib
                alignment_data = AlignIO.read(pathed_aln_name, "fasta") # parse the alignment file
                        
                alignment_sequence_list = list()
                for record in alignment_data:
                    alignment_sequence_list.append(str(record.seq)) # extract just the individual aligned sequences as strings                     

                alignment_matrix = lm.alignment_to_matrix(sequences=alignment_sequence_list, to_type='counts', characters_to_ignore='.-X')
                
            # short alignments fit nicely on a single line
                if alignment_data.get_alignment_length() <= 60:
                    aln_logo = lm.Logo(alignment_matrix, color_scheme='chemistry', figsize=(12,2))
                    aln_logo.ax.set_ylabel('Counts')
                    aln_logo.style_spines(visible=False)
                    aln_logo.style_spines(spines=['left','bottom'], visible=True, linewidth=1)

            # logos with more than 60 stacks are unreadable, so need to be split into multiple rows
                else: 

                # make a temporary full-length logo to determine the appropriate y-axis limits
                    temp_logo = lm.Logo(alignment_matrix, color_scheme='chemistry', figsize=(12,2))
                    ylimits = temp_logo.ax.get_ylim()
                    plt.close()

                # setting row height and alignment slicing start points
                    height_per_row = 2
                            
                    start_slice = 0
                    end_slice = 59
                    counter = 0
                            
                # for alignments whose length is exactly divisible by 60
                    if alignment_data.get_alignment_length() % 60 == 0:
                        rows = alignment_data.get_alignment_length() // 60
                        total_height = height_per_row * rows
                                
                        fig = plt.figure(figsize=(12, total_height), constrained_layout=True)
                        spec = gridspec.GridSpec(ncols=1, nrows=rows, figure=fig)

                        for _ in range(rows):
                            sub_ax = fig.add_subplot(spec[counter, :])
                            aln_logo = lm.Logo(alignment_matrix.loc[start_slice:end_slice], ax=sub_ax, color_scheme='chemistry', figsize=(12,2))
                            aln_logo.ax.set_ylim(ylimits)
                            aln_logo.ax.set_ylabel('Counts')
                            aln_logo.ax.set_anchor('W') # left-justify the plots
                            aln_logo.style_spines(visible=False)
                            aln_logo.style_spines(spines=['left','bottom'], visible=True, linewidth=1)
                            aln_logo.style_xticks(anchor=start_slice, spacing=10)

                            start_slice = end_slice + 1
                            end_slice = end_slice + 60
                            counter = counter + 1
                                    
                # for alignments whose length is not exactly divisible by 60
                    else:
                        rows = (alignment_data.get_alignment_length() // 60) + 1   
                        ratio = (alignment_data.get_alignment_length() % 60) / 60
                        column1 = 12 * ratio
                        column2 = 12 - column1
                        total_height = height_per_row * rows

                        fig = plt.figure(figsize=(12, total_height), constrained_layout=True)
                        spec = gridspec.GridSpec(ncols=2, nrows=rows, width_ratios=[column1, column2], figure=fig)

                        for _ in range(rows-1):
                            sub_ax = fig.add_subplot(spec[counter, :])
                            aln_logo = lm.Logo(alignment_matrix.loc[start_slice:end_slice], ax=sub_ax, color_scheme='chemistry', figsize=(12,2))
                            aln_logo.ax.set_ylim(ylimits)
                            aln_logo.ax.set_ylabel('Counts')
                            aln_logo.ax.set_anchor('W') # left-justify the plots
                            aln_logo.style_spines(visible=False)
                            aln_logo.style_spines(spines=['left','bottom'], visible=True, linewidth=1)
                            aln_logo.style_xticks(anchor=start_slice, spacing=10)

                            start_slice = end_slice + 1
                            end_slice = end_slice + 60
                            counter = counter + 1

                    # draw the short row
                        sub_ax = fig.add_subplot(spec[-1,0])
                        aln_logo = lm.Logo(alignment_matrix.loc[start_slice:], ax=sub_ax, color_scheme='chemistry', figsize=(column1,2))
                        aln_logo.ax.set_ylim(ylimits)
                        aln_logo.ax.set_ylabel('Counts')
                        aln_logo.ax.set_anchor('W') # left-justify the plots
                        aln_logo.style_spines(visible=False)
                        aln_logo.style_spines(spines=['left','bottom'], visible=True, linewidth=1)
                        aln_logo.style_xticks(anchor=start_slice, spacing=10)

            # save logo
                logo_name = base_name + "_logo.png"
                pathed_logo_name = os.path.join(save_path, logo_name)
                plt.savefig(pathed_logo_name, dpi=300)
                plt.close()
                        
                print("Sequence logo saved as %s\n" % logo_name)

                step = step + 1
                break

            else:
                continue

# Step 7: Generate Homology Table
# inputs = 
#   quit, skip step, or continue?
#   FASTA formatted gene list file
# outputs = 
#   * warning if > 100 sequences (require confirmation to continue)
#   tab-delimited homology & taxonomy table file

    while step == 7:

    # check length of input file and throw error message if necessary
        hit_count = 0
        for seq_record in SeqIO.parse(pathed_fastafilename, "fasta"):
            hit_count = hit_count + 1
        if hit_count >= 100:
            print("\nWARNING: Input file contains %s sequences. NCBI prefers that" % str(hit_count))
            print("large requests be performed on weekends or at off-peak times.\n")

    # request clarification about what to do next
        asking = True
        while asking:
            clarify = input("(G)enerate homology and taxonomy table or (S)kip step? ")
            
            if clarify.upper() == "S" or clarify.upper() == "SKIP":
                step = step + 1
                break
            elif clarify.upper() == "Q" or clarify.upper() == "QUIT":
                sys.exit()
            elif clarify.upper() == "H" or clarify.upper() == "HELP":
                helper("homology table")
                continue

            elif clarify.upper() == "G" or clarify.upper() == "GENERATE":
                    
            # generate temporary query file with first record
                first_record = next(SeqIO.parse(pathed_fastafilename, "fasta"))
                SeqIO.write(first_record, "temp.fasta", "fasta")    
    
            # run local bl2seq comparing each hit to the first hit (in that file we just generated)
                xml_filename = name + "_" + datestamp + ".xml"

                print("Running local BLAST comparisons...")
                program = standalone_path + "/blastp"    

                try:
                    cline = NcbiblastpCommandline(cmd=program, query="temp.fasta", subject=pathed_fastafilename, evalue=1, outfmt=5, out=xml_filename) 
                    run = cline()
                except:
                    print("Error! blastp exectutable not found in correct directory.")
                    sys.exit()
    
             # extract useful data and generate table
                print("Generating homology table (requesting taxonomic data from NCBI)...")
                tablefile = name + "_homology_table_" + datestamp + ".txt" 
                pathed_tablefile = os.path.join(save_path, tablefile) 
                save_file = open(pathed_tablefile, "w") 
                first_row = "Species \t Accession \t Length of Homology \t Percent Identical \t Percent Similar \t E Value \t Taxonomy \n"
                save_file.write(first_row)    

                search_handle = open(xml_filename, "r")
                blast_record = SearchIO.read(search_handle, "blast-xml")
                counter = 0
                for hit in blast_record:
                    for hsp in hit:
                        counter = counter + 1
                        species = get_species(hsp.hit_description)            
                        accession = hsp.hit_id              
                        ident = (hsp.ident_num * 100) / hit.seq_len
                        sim = (hsp.pos_num * 100) / hit.seq_len
                        e = hsp.evalue

                    # generate taxonomy column, which has historically been the most likely part of this script to crash  
                        handle = Entrez.esearch(db = "taxonomy", term = species)
                        record = Entrez.read(handle)
                        try:
                            tax_id = record["IdList"][0]
                        except:
                            tax_id = "Not available"
                            
                        if tax_id != "Not available":
                            handle = Entrez.efetch(db = "taxonomy", id = tax_id)
                            records = Entrez.read(handle)
                            try:
                                taxon = records[0]["Lineage"]
                            except:
                                taxon = "Not available"
                        else:
                            taxon = "Not available"
                            
                    # write the row of extracted data into the table
                        print("Table row: " + str(counter))
                        row = species + "\t" + accession + "\t" + str(hit.seq_len) + " / " + str(len(first_record)) + "\t" + str(ident) + "\t" + str(sim) + "\t" + str(e) + "\t" + taxon + "\n"
                        save_file.write(row)

                print("Homology table of data for " + str(counter) + " aligned regions saved as " + tablefile + "\n")
                save_file.close()
                search_handle.close()

                # eliminate temporary file and XML file
                os.remove('temp.fasta')
                os.remove(xml_filename)

                asking = False
                step = step + 1
                break

            else:
                continue

# Give the option to end the program, to repeat an analysis, or to start over with a new protein
    while step >= 8:

    # request clarification about what to do next
        asking = True
        while asking:
            clarify = input("(R)epeat analysis, analyze a (N)ew protein, or (Q)uit? ") 
            if clarify.upper() == "R" or clarify.upper() == "REPEAT" or clarify.upper() == "REPEAT ANALYSIS":
                step = 2
                break
            elif clarify.upper() == "H" or clarify.upper() == "HELP":
                helper("repeat step")
                continue
            elif clarify.upper() == "N" or clarify.upper() == "NEW" or clarify.upper() == "NEW PROTEIN":
                step = 1
                break
            elif clarify.upper() == "Q" or clarify.upper() == "QUIT":
                sys.exit()
            else:
                continue
