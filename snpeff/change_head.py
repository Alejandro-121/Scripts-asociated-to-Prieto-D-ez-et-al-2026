import os

os.chdir("/home/alejandro/Documentos/nextcloud/samoa/data/S288C_reference_genome_R64-1-1_20110203")

in_fasta = "S288C_reference_sequence_R64-1-1_20110203.fsa"
out_fasta = "new_head.fasta"


file = open(in_fasta, "r")
fasta = file.readlines()


# new headers
headers = ("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII",
           "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV",
           "chrXV", "chrXVI", "chrMito")



# change headers
ite = 0
for i in range(len(fasta)):
    if fasta[i][0] == ">":
        print("old head = " + fasta[i])
        fasta[i] = ">" + headers[ite]
        ite += 1

# CAN SOMEONE EXPLAIN WHY THIS DOES NOT WORK!!!!!!!!!!
#ite = 0
#for line in fasta:
#    if line[0] == ">":
#        print("old head = " + line)
#        line = ">" + headers[ite]
#        ite += 1


# check new head
print("NEW HEAD STARTS HERE!!!!")
for i in range(len(fasta)):
    if fasta[i][0] == ">":
        print("NEW head = " + fasta[i])


# write output
with open(out_fasta, "w") as file:

    for line in fasta:
        if line[0] == ">":
            file.write(line + "\n")
        else:
            file.write(line)
print("File written")
# Reference genome chosen
S288C_reference_genome_R64-1-1_20110203.tgz
