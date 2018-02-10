#### Import libraries ####
import math, time

#### Import files #####
genbank_fn = "homework5_example_2.gbff"

#### Classes ####
class CDS(object):
    def __init__(self,raw):
        self.raw = raw
        #set strand
        if raw[0] == 'c': #complement will always come first
            self.strand = 1 #reverse
        else:
            self.strand = 0 #forward

        self.coords = [[int(inner) for inner in x.split('..')] for x in raw.replace('(',',').replace(')',',').split(',') if x not in ("","join","complement")]
    def get21(self):

        if self.strand == 0: #forward
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
                tss21 = seq[self.coords[-2][1]-len_left:self.coords[-2][1]+1] + tss21_end# theoretically could hit something where there is no exon 2 and a super short e1...
            return reverse_complement(tss21)
    def pos(self):
        if self.strand == 0: #forward
            start = self.coords[0][0]
        else:
            start = self.coords[-1][1]
        return start

    def score(self):
        the_tss21 = self.get21()
        score = sum( weight_matrix[i][the_tss21[i]] for i in range(21)) #sum across the weights
        return score

#for getting a "CDS" from the position
class CDS_pos(object):
    def __init__(self,pos,strand):
        self.position = pos
        self.strand = strand #0 = forward, 1 = reverse

    def get21(self):
        start = self.position
        tss21 = seq[start-10:start+11]
        if self.strand == 0: #forward
            return tss21
        else: #reverse
            return reverse_complement(tss21)

    def score(self):
        the_tss21 = self.get21()
        score = sum(weight_matrix[i][the_tss21[i]] for i in range(21)) #sum across the weights
        return score


#### Funcions ####
def reverse_complement(sequence): #reverse complement
    comp_dict = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return ''.join(comp_dict[base] for base in sequence[::-1]) #reverse and complement

def read_gb(filename):
    try:
        f = open(filename,'r')
    except IOError:
        print("The file %s does not exist" % filename)
        return
    else:
        file_list = [line.strip() for line in f if line.strip().split()[0] != "/"]
        more_lines = False
        cds_multline = []
        holder = ""
        for line in file_list: #Man, I really hate genbank format...
            if line.split()[0] == 'CDS': #first line of CDS
                holder = line.rstrip()
                more_lines = True
            elif '=' in line and more_lines == True: #line after CDS
                cds_multline.append(holder)
                more_lines = False
                holder = ""
            elif ")" not in line and more_lines == True: #middle line(s) of CDS
                holder += line.rstrip()
            elif ")" in line and more_lines == True: #last line of CDS
                holder += line.rstrip()
                cds_multline.append(holder)
                more_lines = False
                holder = ""

        cds_list = [line.strip().split()[1] for line in cds_multline if "<" not in line and ">" not in line] #line.strip().split()[0] == 'CDS' and
        #ADDED 'n' TO BEGINNING of sequence so the indices will now match the genbank ones :D
        sequence = "n" + ''.join([''.join(line.strip().split()[1:]) for line in file_list if line.strip().split()[0].isdigit()])
        return sequence.lower(), cds_list

def getKey(item): #for sorting things based on first
    return int(item[0])

#### Code ####
t = time.time()

seq, raw_cds_list = read_gb(genbank_fn) #parse seq and unprocessed cds lines

cds_list = [CDS(raw) for raw in raw_cds_list] #put the cds into an object

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

wrong_cds = []
wrong_start = []

cds_starts = []

for gene in cds_list:
    tss21 = gene.get21()
    cds_starts.append((gene.pos(),gene.strand))
    if tss21[10:13] != "atg":
        print(tss21)
        wrong_cds.append(gene.raw) #double check I got the right sequences; should all be atg
        wrong_start.append(tss21[10:13])

    for i in range(21):
        basecounts[i][tss21[i]] += 1

counter = -10
for i in range(21):
    print(counter, basecounts[i]["a"],basecounts[i]["c"],basecounts[i]["g"],basecounts[i]["t"])
    counter += 1

basecounts_no_n = [{nt: basecounts[i][nt] for nt in ['a','c','g','t']} for i in range(21)]

basefreqs = [ {nt: basecounts_no_n[i][nt]/sum(basecounts_no_n[i].values()) for nt in ['a','c','g','t']} for i in range(21)]

# Print frequency matrix
print("")
print("Frequency Matrix:")
counter = -10
for i in range(21):
    print(counter,
          format(basefreqs[i]["a"],'.4f'),
          format(basefreqs[i]["c"],'.4f'),
          format(basefreqs[i]["g"],'.4f'),
          format(basefreqs[i]["t"],'.4f')
          )
    counter += 1

# Print Weight Matrix
print("")
print("Weight Matrix:")

weight_matrix = []
for i in range(21):
    holder_dict = {}
    for nt in ['a','c','g','t','n']:
        if nt == 'n':
            holder_dict[nt] = 0
        elif basefreqs[i][nt] != 0:
            holder_dict[nt] = math.log2( basefreqs[i][nt] / background_freq[nt] )
        else:
            holder_dict[nt] = -99
    weight_matrix.append(holder_dict)

counter = -10
for i in range(21):
    print(counter,
          format(weight_matrix[i]["a"],'.4f'),
          format(weight_matrix[i]["c"],'.4f'),
          format(weight_matrix[i]["g"],'.4f'),
          format(weight_matrix[i]["t"],'.4f')
          )
    counter += 1

#Print Maximum Score
print("")
max_score = sum(max(weight_matrix[i].values()) for i in range(21))
print("Maximum Score:",format(max_score,'.10f'))

#Print Score Histogram CDS
print("")
print("Score Histogram CDS:")

#cds_score_set = {math.floor(gene.score()) for gene in cds_list} #set comprehension

cds_score_list = [math.floor(gene.score()) for gene in cds_list]
cds_score_hist = [(score, cds_score_list.count(score)) for score in sorted(set(cds_score_list))]
for score_line in cds_score_hist:
    print(score_line[0],score_line[1])


#Print Score Histogram All
print("")
print("Score Histogram All:")

all_sites = [CDS_pos(x,strand) for strand in [0,1] for x in range(11,len(seq) - 10)]
all_sites_21s = {(gene.position,gene.strand) : gene.score() for gene in all_sites}

all_score_list = [math.floor(val) if math.floor(val) >= -50 else -50 for val in all_sites_21s.values()]
all_score_set = sorted(list(set(all_score_list))) #list of sorted set
all_score_hist = [(score, all_score_list.count(score)) for score in all_score_set[1:]]
for score_line in all_score_hist:
    print(score_line[0],score_line[1])
print("lt-50",all_score_list.count(-50))

# Print Position List
print("")
print("Position List:")

#print(len(all_sites_21s), len(cds_starts))

post_list = []
for key, value in all_sites_21s.items():
    if value > 10 and key not in cds_starts:
        post_list.append((key[0],key[1], value))

for (a,b,c) in sorted(post_list, key=getKey):
    print(a,b,format(c,'.4f'))




print(time.time() - t, " seconds")