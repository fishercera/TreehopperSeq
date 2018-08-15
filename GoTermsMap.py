# ~~~~ Cera Fisher (2018)
# Takes a file of gene IDs and GoTerms from EnTap and makes a many-to-many tab delimited mapping file appropriate for use with GoSeq in R.
# Create the infile in bash using:
# cut -f 1,28,29,30 entap.final_annotations.tsv > entap_goTerms.tsv
import sys

infile=sys.argv[1]
outfile=sys.argv[2]

# This function takes as input the gene id, one of the three GOTerm fields from EnTAP, and an output file
# (specified as an argument when the script is run)
# It splits the list of GoTerms on the characters "GO:"
# (rather than the comma delimiter, because some GOTerms have commas in them)
# and writes one line with the gene id for each goterm

def writeGoLine(id, field, outfile):
    if field != '':                                     # We don't want to bother with blank fields.
        Terms = field.split("GO:")                      # Splitting the field into an array, at the character "GO:" bc some GOterms have commas
        for term in Terms:                              # Going to loop through this newly created array
            term = term.strip()                         # strip spaces off both sides
            term = term.rstrip(",")                     # strip the trailing comma
            if term != '':                              # use this to skip the incidentally created blank first item
                with open(outfile, "a") as mapping:     # writing with with open bc it's better for memory
                    mapping.write(id+'\tGO:'+term+'\n') # Sneakily adding back in the "GO:" that we lost when we split the field on it


with open(infile, "r") as goterms: # Open the infile to read
    for line in goterms:                                # This is where we read the file created from EnTAP output
        Cols = line.split("\t")                         # Split apart the different fields
        id = Cols[0].strip()                            # Get the transcript ID
        goBio = Cols[1].strip()                         # Get the list of Biological Process terms
        goCell = Cols[2].strip()                        # And the list of Cellular Component Terms
        goMol = Cols[3].strip()                         # And the list of Molecular function terms
                                        # Below, if you want ONLY one kind of annotation (bio process, cell. component)
                                        # comment out the lines that write the kind of annotation you DON'T want.
        writeGoLine(id, goBio, outfile)                 # Here we use the function defined above to write the one-to-one mapping file.
        writeGoLine(id, goCell, outfile)                # Note, you could modify this to specify three different output files, and 
        writeGoLine(id, goMol, outfile)                 # write each type of GO term to a separate file, but it's very unnecessary. 


print("I think I'm done")                               # I like when my scripts talk to me. 

# IMPORTANT NOTE: the output file may still need to be modified in order to use it with goseq. GoSeq wants only the GO identifier, e.g.
# GO:0032991 but our EnTAP annotation is a bit more verbose, listing the GO identifier and the term, a la
#               GO:0003008-system process(L=2)

# For my own purposes, I decided to take the output file made by this script and, in bash, do
# cut -f 2 outputfile.txt (cutting the GOterm column) > fileofjustgoterms.txt
# sed -i "s/-.*$//g" (using regular expressions to find and replace the text that occurs right after the GO id, which starts with a hyphen
# paste outputfile.txt fileofjustgorterms.txt > newthreecolumnfile.txt
