#### Import files #####
genbank_fn = "homework5_example_2.gbff"

#### Classes ####
class CDS(object):
    def __init__(self,raw):
        self.raw = raw
        #set strand
        if raw[0] == 'c': #complement will always come first
            self.strand = "reverse"
        else:
            self.strand = "forward"

        self.coords = [[int(inner) for inner in x.split('..')] for x in raw.replace('(',',').replace(')',',').split(',') if x not in ("","join","complement")]
    def get21(self):

        if self.strand == "forward":
            start = self.coords[0][0]
            exon1_len = self.coords[0][1] - start
            if exon1_len >= 10 and start > 10 and start + 10 < len(seq): #if the right length and
                tss21 = seq[start-10:start+11] #slice is exclusive at end
            elif start < 10 or start + 10 >= len(seq): #too short for matrix
                tss21 = "n" * 21 #"throw out" sequences that are too short
            else: #
                try:
                    tss21 = seq[start-10:self.coords[0][1]+1] #first bit
                    len_left = 21 - len(tss21)
                    tss21 += seq[self.coords[1][0]:self.coords[1][0] + len_left] #theoretically could hit something where there is no exon 2 and a super short e1...
                except:
                    print(self.coords, exon1_len)
                    raise
            #the part in exon "2"
            return tss21
        else: #reverse
            start = self.coords[-1][1]
            #print("start", start)
            exon1_len = start - self.coords[-1][0]
            if exon1_len >= 10 and start > 10 and start + 10 < len(seq):  # SAME
                tss21 = seq[start - 10:start + 11]  # slice is exclusive at end
            elif start < 10 or start + 10 >= len(seq):  # SAME
                tss21 = "n" * 21  # "throw out" sequences that are too short
            else:  #
                tss21_end = seq[self.coords[-1][0]: start + 11]
                len_left = 21 - len(tss21_end)
                tss21 = seq[self.coords[-2][1]-len_left:self.coords[-2][1]] + tss21_end# theoretically could hit something where there is no exon 2 and a super short e1...
            return reverse_complement(tss21)


#### Funcions ####
def reverse_complement(sequence): #reverse complement
    comp_dict = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return ''.join(comp_dict[base] for base in sequence[::-1]) #reverse and complement

def read_gb(filename):
    try:
        f = open(filename,'r')
    except IOError:
        print("The file %s does not exist" % filename)
        return
    else:
        file_list = [line.strip() for line in f if line.strip().split()[0] == "CDS" or line.strip().split()[0].isdigit()]
        cds_list = [line.strip().split()[1] for line in file_list if line.strip().split()[0] == 'CDS' and "<" not in line and ">" not in line]
        #ADDED 'n' TO BEGINNING of sequence so the indices will now match the genbank ones :D
        sequence = "n" + ''.join([''.join(line.strip().split()[1:]) for line in file_list if line.strip().split()[0].isdigit()])
        return sequence.lower(), cds_list

seq, raw_cds_list = read_gb(genbank_fn) #parse seq and unprocessed cds lines

cds_list = [CDS(raw) for raw in raw_cds_list] #put the cds into an object

#test = CDS("complement(join(1..11,15..19))")
#test = CDS("join(12..17,22..29)")

#print(cds_list[0].raw)
#print(cds_list[0].get21())

#Print header
print("Assignment: GS540 HW5")
print("Name: Emma Hoppe")
print("Email: erhoppe@uw.edu")
print("Language: Python3")
print("Runtime: ")
print("")

# Print the Nucleotide Histogram
total_count_dict = {}
for nt in seq[1:]: #since I added the n to the beginning to get the indices right
    total_count_dict[nt] = total_count_dict.get(nt,0)+1
print("Nucleotide Histogram:")
for base in ["a","c","g","n","t"]:
    print(base.upper(),"=",total_count_dict[base])

#Print Background Frequency
# We need to grab the
print("")
print("Background Frequency:")
total_non_ambig = sum(total_count_dict.values()) - total_count_dict["n"]
background_freq = {
    "a": (total_count_dict["a"] + total_count_dict["t"]) / (2*total_non_ambig),
    "c": (total_count_dict["c"] + total_count_dict["g"]) / (2*total_non_ambig),
    "g": (total_count_dict["g"] + total_count_dict["c"]) / (2*total_non_ambig),
    "t": (total_count_dict["t"] + total_count_dict["a"]) / (2*total_non_ambig),
}
for base in ["a","c","g","t"]:
    print(base.upper(),"=",round(background_freq[base],4))

#Print Count Matrix
print("")
print("Count Matrix:")
basecounts = [{'a': 0, 'c': 0, 'g': 0, 't':0, 'n': 0} for i in range(21)] #initialize the dictionary
for gene in cds_list:
    tss21 = gene.get21()
    if tss21[10:13] is not "atg":
        print(gene.raw)

    for i in range(21):
        basecounts[i][tss21[i]] += 1

counter = -10
for i in range(21):
    print(counter, basecounts[i]["a"],basecounts[i]["c"],basecounts[i]["g"],basecounts[i]["t"],)
    counter += 1

### some of these are off; should print all of the first three letters which will all be atg

