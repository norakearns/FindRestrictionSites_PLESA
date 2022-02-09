# Packages
from Bio import SeqIO
import re
import itertools
from itertools import product

# Files
cuts_file = open("RE_cutsites.txt","wt")
cut_sites = open("cut_sites.txt", "r") # file containing Enzymes and their known cut sites
oligos = open("oligos.fasta", "rt")
Cut_locations = open("cut_locations.txt","wt")


'''
OBJECTIVE:
Figure out all AA combinations that can make a restriction site. Locate those cut sites in an input AA sequence.
Input: AA sequence
Output: potential RE cut sites and their positions

METHOD:
 read in Enzyme cut site:
    make frame 1: sequence as is (GAATTC)
        split the string into 123 | 456 (GAA | TTC)
        use dictionary to get 1-letter AAs (E | F)
    
    make frame 2: N + sequence + NN (NGAATTCNN)
        split the string into 123 | 456 | 789 (NGA | GTT | CNN)
        make 4 strings, subbing N for every possible nucleotide (ex: UGA, TGA, CGA, AGA)
        use dictionary to get 1-etter AAs
        make all possible two letter combinations and generate 16 codon possibilities 
        translate codon possibilities to 1-letter AAs
        make all possible AA strings
    
    make frame 3: NN + sequence + N (NNGAATTCN)
    
'''

# Dictionary of all codons and their one-letter Amino Acid translations
AA_dict = {'TTT':'F', 'TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*', 'TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
           'CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R','ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K',
           'AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

# Example sequence
Seq = "ETVVLVHRVMAVFNAYLLVFTIIWLRARLLQPWRQLLAMASAVSHRDFTQRANISGRNEMAMLGTALNNMSAELAESYAVLEQRVQEKTAGLEHKNQILSFLWQANRRLHS"

def reverse_translate(AA_seq):
 '''
 create all possible reverse translated DNA sequences from an amino acid string
 input: NAY
 output: ['AATGCTTAT', 'AATGCTTAC', 'AATGCCTAT', 'AATGCCTAC', 'AATGCATAT', 'AATGCATAC', 'AATGCGTAT', 'AATGCGTAC', 'AACGCTTAT', 'AACGCTTAC', 'AACGCCTAT', 'AACGCCTAC', 'AACGCATAT', 'AACGCATAC', 'AACGCGTAT', 'AACGCGTAC']
 '''
    possible_codons = []
    for i in AA_seq:
        codons = [k for k, v in AA_dict.items() if v == i]
        possible_codons.append(codons)
    possible_strings_iters = list(itertools.product(*possible_codons))
    possible_strings = []
    for i in range(0, len((possible_strings_iters))):
        possible_strings.append(''.join(possible_strings_iters[i]))
    return possible_strings

set = ['R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N']
def containsAny(str, set):
    '''
    Find if any weird characters in the enzyme cut site
    :param str: ACRT
    :param set: ['R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N']
    :return: 1 (0 if R not in string)
    '''
    for c in set:
        if c in str: return 1
    return 0

sites_dict = {}
for line in cut_sites:
# Create a dictionary of the enzyme cut sites (sites_dict)
    enzyme = cut_sites.readline().split('	')[0] # read in the cut_sites file and grab the enzyme name
    site = cut_sites.readline().split('	')[1].strip('\n') # read in the cut_sites file and grab the sequence
    if containsAny(site, ['R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D','r','y','m','k','s','w','h','b','v','d']) == 0: # if a sequence has any weird characters
        sites_dict[enzyme]=site.upper() # store only the normal sequences
    else:
        continue

# Dictionary of enzymes and their cut sites
sites_dict = {'AarI': 'GACNNNNNNGTC', 'Aba6411II': 'GCAAAC', 'AccI': 'CGCG',
              'Acc36I': 'GGTACC', 'AccB1I': 'CCANNNNNTGG', 'AceIII': 'AGCCAG',
              'AclI': 'GGATC', 'AcuI': 'CACGTG', 'AdeI': 'GAANCAG', 'AfeI': 'CCNNNNNNNGG',
              'AflIII': 'ACCGGT', 'AhaIII': 'GACNNNNNGTC', 'AjnI': 'GAANNNNNNNTTGG', 'AlfI': 'GAACNNNNNNTCC', 'AluBI': 'GGATC', 'Alw26I': 'GTGCAC', 'Aod1I': 'TCCGGA', 'AoxI': 'GGGCCC', 'ApyPI': 'GCCGNAC', 'AseI': 'GATC', 'Asp337I': 'GAANNNNTTC', 'AspA2I': 'ACCCAC', 'AspDUT2V': 'CGCCCAG', 'AspNIH4III': 'GGNCC', 'Asp114pII': 'GGNCC', 'AvaII': 'ATGCAT', 'Awo1030IV': 'CCTNAGG', 'BaeGI': 'CCCGAG', 'BanLI': 'GAAGNNNNNNTAC', 'Bau1417V': 'GGCGAG', 'BbvI': 'GAAGAC', 'BbvCI': 'CCATC', 'Bce3081I': 'TATCNAG', 'BceSIV': 'ACGGC', 'BciT130I': 'GTATCC', 'BcuI': 'TGANNNNNNTCA', 'BfaI': 'GANGGAG', 'BfuI': 'ACCTGC', 'BglI': 'AGATCT', 'BlnI': 'GAGGAC', 'BlsI': 'AGTACT', 'Bme1390I': 'C', 'BmgT120I': 'GGNNCC', 'BmrFI': 'GCATC', 'BmuI': 'GACNNNNGTC', 'BplI': 'CTGGAG', 'Bpu10I': 'TTCGAA', 'BsaHI': 'CCNNGG', 'BsaXI': 'CAACAC', 'BscAI': 'CCCGT', 'Bse1I': 'GATNNNNATC', 'Bse118I': 'TCCGGA', 'BseCI': 'CCNNGG', 'BseGI': 'GATNNNNATC', 'BseMI': 'CTCAG', 'BsePI': 'GAGGAG', 'BseXI': 'CGGCCG', 'BsgI': 'CGCG', 'BshVI': 'CACGAG', 'BsiWI': 'CCNNNNNNNGG', 'BslFI': 'GAATGC', 'BsmBI': 'GGGAC', 'Bsp19I': 'GACNNNNNNTGG', 'Bsp119I': 'GGGCCC', 'Bsp1720I': 'CCGCAT', 'BspANI': 'CTCAG', 'BspD6I': 'TCCGGA', 'BspGI': 'TCATGA', 'BspLU11I': 'ACCTGC', 'BspMAI': 'CCAGA', 'BspPI': 'GCTCTTC', 'BsrI': 'CCGCTC', 'BsrFI': 'TGTACA', 'BssECI': 'GCGCGC', 'BssNI': 'GTATAC', 'BssT1I': 'CTCTTC', 'BstACI': 'CTTAAG', 'BstAUI': 'TTCGAA', 'BstBAI': 'ACNGT', 'BstENI': 'GGATG', 'BstH2I': 'GCGC', 'BstMAI': 'GATC', 'BstPI': 'GACNNNNGTC', 'BstV2I': 'CCANNNNNNTGG', 'BstYI': 'CGGCCG', 'BsuI': 'ATCGAT', 'BsuRI': 'ATCGAT', 'BtgZI': 'GCNGC', 'BtsI': 'CAGTG', 'BtuMI': 'ACCTGC', 'Cac8I': 'CAGNNNCTG', 'CchII': 'CCCAAG', 'CciI': 'GCGGCCGC', 'CcrNAIII': 'CATCG', 'Cdi13746V': 'CCGATCC', 'CfrI': 'CCCGGG', 'Cfr13I': 'CCGCGG', 'Cgl13032II': 'GATC', 'Cko11077IV': 'ATCGAT', 'Csa9238II': 'GACGC', 'CspI': 'GTAC', 'CspAI': 'CAANNNNNGTGG', 'CviQI': 'TGCA', 'Dde51507I': 'GGCGCC', 'DpnII': 'TTTAAA', 'DraIII': 'CAAGNAC', 'DrdII': 'TACGAC', 'DsaI': 'GACNNNNNNGTC', 'Eam1104I': 'GACNNNNNGTC', 'EciI': 'GAGCTC', 'Ecl35734I': 'CGGCCG', 'Eco31I': 'GATATC', 'Eco47III': 'CGGCCG', 'Eco72I': 'CCTNAGG', 'Eco91I': 'TACGTA', 'Eco147I': 'GCACAG', 'Eco8164I': 'GAAANTC', 'EcoHSI': 'GAGCTC', 'EcoMVII': 'CCTNNNNNAGG', 'EcoRII': 'GATATC', 'EgeI': 'GGCGCC', 'EsaBC3I': 'GACCAC', 'Esp3I': 'CAGAAG', 'FaiI': 'AAGNNNNNCTT', 'FatI': 'CCCGC', 'FbaI': 'AGAAGG', 'Fco1691IV': 'GGGAC', 'FnuDII': 'GCNGC', 'FriOI': 'GGCCGGCC', 'FspAI': 'CTAG', 'GauT27I': 'ATGCAC', 'GlaI': 'GCNGC', 'GsaI': 'CTGGAG', 'HaeII': 'GGCC', 'HapII': 'TGGCCA', 'Hca13221V': 'CGANNNNNNTCC', 'Hin1I': 'CATG', 'Hin4II': 'GCGC', 'HinfI': 'GTTAAC', 'HphI': 'GTNNAC', 'Hpy99XXII': 'GTNNAC', 'Hpy188I': 'TCNNGA', 'HpyAV': 'GCGTA', 'HpyAXVI-mut2': 'GGANNAG', 'HpyCH4III': 'ACGT', 'HpyF3I': 'GCNNNNNNNGC', 'HpyPU007XIX': 'ACGT', 'HpyUM032XIII-mut1': 'GAAAG', 'HspAI': 'GAGCAGC', 'KpnI': 'TCCGGA', 'Kpn327I': 'GNGCGAG', 'KpnNH25III': 'GTTCNAC', 'KspI': 'TGATCA', 'KspAI': 'GATC', 'Lbr124II': 'ACAAAG', 'LlaG50I': 'GCTCC', 'Lpn11417II': 'GCNCAAC', 'Lra68I': 'TGGAAT', 'MabI': 'CTAG', 'MaeIII': 'GATC', 'MauBI': 'AGGCGA', 'MboI': 'GAAGA', 'Mcr10I': 'CAATTG', 'MhlI': 'GTNNAC', 'Mla10359I': 'TGGCCA', 'Mlu211III': 'AATT', 'MlyI': 'GGCGCC', 'MnlI': 'TGGCCA', 'MreI': 'TCCGGA', 'MroXI': 'TGGCCA', 'MslI': 'CCGG', 'MspA1I': 'CTTAAG', 'MspJI': 'CCNGG', 'MssI': 'TGCGCA', 'MtuHN878II': 'CAATTG', 'Mva1269I': 'CGCG', 'NaeI': 'ACCAGC', 'NdeI': 'GATC', 'NgoAVIII': 'GCCGGC', 'NlaIV': 'CATCAC', 'NmeAIII': 'GCCGAC', 'NmuCI': 'GCGGCCGC', 'NruI': 'TGCGCA', 'NspI': 'TTCGAA', 'NspES21II': 'ACGAG', 'PacI': 'GTAATC', 'PaePA99III': 'CTCGAG', 'Pal408I': 'GGCGCGCC', 'PasI': 'GCGCGC', 'Pbu13063II': 'GACGAG', 'PciI': 'GCTCTTC', 'PcsI': 'GAATGC', 'Pdi8503III': 'GAANNNNTTC', 'Pfl1108I': 'CCCTNAG', 'Pfl10783II': 'GACNNNGTC', 'PflPt14I': 'TCCNGGA', 'PfrJS12V': 'CTTCNAC', 'PkrI': 'CATCAG', 'Ple19I': 'CGCCGAC', 'PmaCI': 'GTTTAAAC', 'PmlI': 'GAACNNNNNCTC', 'PpsI': 'ATGCAT', 'PpuMI': 'CAGANGC', 'Pse18267I': 'GACNNNNGTC', 'Psp6I': 'GCGAAG', 'Psp124BI': 'CACGTG', 'PspEI': 'CCCAGC', 'PspLI': 'GGNNCC', 'PspOMII': 'GGNCC', 'Pst145I': 'GATCGAG', 'PsyI': 'GCGCGC', 'PvuII': 'GAAAGAG', 'RceI': 'CCGCAG', 'RlaII': 'CCCACA', 'RpaI': 'CCCGCAG', 'RpaTI': 'TCGCGA', 'Rsp008V': 'CACACG', 'SacI': 'CCGCGG', 'Sag901I': 'GTCGAC', 'SapI': 'TTAA', 'SauI': 'GGNCC', 'Sau5656II': 'GTANNNNNNTGG', 'Sbo46I': 'AGTACT', 'SciI': 'GCTAAT', 'Sdy9603I': 'CCNNGG', 'Sen5794III': 'GTTCAT', 'SfaAI': 'GCATC', 'SfeI': 'GGCCNNNNNGGCC', 'SfoI': 'CTCGAG', 'SfuI': 'CNNG', 'SgrBI': 'CGTCGACG', 'SgsI': 'GGGTC', 'SlaI': 'CCCGGG', 'Sma10259II': 'CTTGAC', 'SmoI': 'GTATAC', 'SnaBI': 'GGCCGAG', 'Spe19205IV': 'GCATGC', 'Sse9I': 'CGCCGGCG', 'Sse8647I': 'AGGCCT', 'SspI': 'CGCAGCG', 'SspDI': 'GGTGA', 'SspMI': 'GAGCTC', 'Sth132I': 'CCGG', 'SthSt3II': 'GGATG', 'StyI': 'CCNGG', 'SwaI': 'ACNGT', 'TaqI': 'GACCGA', 'TfiI': 'GTGAAG', 'TpyTP2I': 'TTAA', 'Tsp4CI': 'ATGAA', 'TspGWI': 'CCCGGG', 'TssI': 'CACNNNNNNTCC', 'UbaF11I': 'CTACNNNGTC', 'UbaF14I': 'CGAACG', 'UnbI': 'CCANNNNNTGG', 'XmaIII': 'CCTAGG', 'XmnI': 'CTAG', 'Yps3606I': 'AGGAAG', 'ZrmI': 'ATGCAT'}

t2 = list(product('ATGC',repeat=2)) # create all possible combinations of 2 NT 
t3 = list(product('ATGC',repeat=3)) # create all possible combinations of 3 NT 
twoNT_list=[]
threeNT_list=[]
for i in range(0,len(t2)):
    twoNT_list.append(''.join(t2[i])) # create all possible 2 NT strings
for i in range(0,len(t3)):
    threeNT_list.append(''.join(t3[i])) # create all possible 3 NT strings

# all possible codon patterns with Ns in them 
end_one_N = re.compile("[ACTG]{2}[N]")
begin_one_N = re.compile("[N][ACTG]{2}")
end_two_N = re.compile('[ACTG][N]{2}')
begin_two_N = re.compile("[N]{2}[ACTG]")
flanking_N = re.compile("[N][ACTG][N]")
middle_N = re.compile("[ACTG][N][ACTG]")
begin_2_N = re.compile("[ACTG][N]")
end_2_N = re.compile("[N][ACTG]")

# C1_F1 = codon 1, frame 1
for enzyme in sites_dict.keys():
    frame_1 = sites_dict[enzyme] # GAATTC
    frame_2 = "N" + sites_dict[enzyme] + "NN" # NNGAATTCN
    frame_3 = "NN" + sites_dict[enzyme] + "N" # NGAATTCNN
    sites = [frame_1, frame_2, frame_3] # [GAATTC, NNGAATTCN, NGAATTCNN]
    print(sites)
    Nucleotides = ["A","C","T","G"]
    for site in sites:
        frame1_codons = [site[i:i + 3] for i in range(0, len(site), 3)] # break the sequence at every 3rd NT
        if len(frame1_codons[-1]) == 1: # if the length of the last codon is 1, add two NN on the end
            frame1_codons[-1] = frame1_codons[-1] + "NN" 
        if len(frame1_codons[-1]) == 2: # if the last codon is only 2 long, add all potential nucleotides to the end to make triplets
            frame1_codons[-1] = frame1_codons[-1] + "N"

        for index in range(len(frame1_codons)):
            if type(frame1_codons[index]) == str and frame1_codons[index]=="NNN":
                frame1_codons[index] = threeNT_list
            elif type(frame1_codons[index]) == str and re.match(begin_two_N, frame1_codons[index]):
                possible_codons_2 = []
                for twoNTs in twoNT_list:
                    possible_codons_2.append(twoNTs + frame1_codons[index][-1])
                frame1_codons[index] = possible_codons_2 # add all possible combinations of 2 NTs to the beginning of the sequence
            elif type(frame1_codons[index]) == str and re.match(end_two_N, frame1_codons[index]) and len(frame1_codons[index]) ==3:
                possible_codons_2 = []
                for twoNTs in twoNT_list:
                    possible_codons_2.append(frame1_codons[index][0] + twoNTs)
                frame1_codons[index] = possible_codons_2
            elif type(frame1_codons[index]) == str and re.match(begin_one_N, frame1_codons[index]) and len(frame1_codons[index])==3:
                possible_codons_2 = []
                for i in Nucleotides:
                    new_string = str(i) + str(frame1_codons[index][-2:])
                    possible_codons_2.append(new_string)
                frame1_codons[index] = possible_codons_2 # add all possible combinations of 2 NTs to the beginning of the sequence
            elif type(frame1_codons[index]) == str and re.match(end_one_N, frame1_codons[index]):
                possible_codons_2 = []
                for NT in Nucleotides:
                    new_string = str(frame1_codons[index][:2]) + str(NT)
                    possible_codons_2.append(new_string)
                frame1_codons[index] = possible_codons_2
            elif type(frame1_codons[index]) == str and re.match(flanking_N, frame1_codons[index]):
                possible_codons_2 = []
                for NT in Nucleotides:
                    new_string_2 = str(NT) + str(frame1_codons[index][-2:])
                    for NT in Nucleotides:
                        new_string_3 = new_string_2[:2] + NT
                    possible_codons_2.append(new_string_3)
                frame1_codons[index] = possible_codons_2
            elif type(frame1_codons[index])==str and re.match(middle_N, frame1_codons[index]):
                possible_codons = []
                for NT in Nucleotides:
                    new_string = str(frame1_codons[index][0]) + NT + str(frame1_codons[index][2])
                    possible_codons.append(new_string)
                frame1_codons[index] = possible_codons_2
            elif type(frame1_codons[index]) == str:
                possible_codons = []
                possible_codons.append(frame1_codons[index])
                frame1_codons[index] = possible_codons

        frame1_AAs=[]
        for codon_array in frame1_codons:
        # translate codons to one letter AAs
            letter_array = []
            for codon in codon_array:
                one_letter_AA = AA_dict[codon]
                if one_letter_AA not in letter_array:
                    letter_array.append(one_letter_AA)
            frame1_AAs.append(letter_array)

        possible_strings_iters = list(product(*frame1_AAs)) # create strings of all possible cut sites
        possible_strings = []
        for i in range(0, len((possible_strings_iters))):
            possible_strings.append(''.join(possible_strings_iters[i]))

        cuts_file.write(enzyme + ":" + site + ":" + str(len(possible_strings)) + "\n")
        for i in possible_strings:
            print(i)
            cuts_file.write(i + "\n")
            site = Seq.find(i) # search for cut sites in a sequence
            if site != -1:
                Cut_locations.write(enzyme + "\n")
                Cut_locations.write(str(site) + ", ")
                Cut_locations.write(i + ", ")
                Cut_locations.write(str(reverse_translate(i)) + "\n")





# i want the enzyme name, position start in a.a., position start frame (a,b,c,), position end in a.a., position end frame, codon combination that allows that cut


