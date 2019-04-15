import sys

infilename = sys.argv[1]
outfilename = sys.argv[2]

taxon_seq_dict = {}
num_heliconius = 0
num_nymphalid = 0
e_and_a = 0

heliconius_list = ["HeraHhimHyb", "Hhim", "HeraDisco", "HeraRef",
                    "Hhsa", "Hdem", "Hsar", "Htel", "HmelRef", "HmelDisco",
                    "Htim", "Hcyd", "Hpar", "Hbes", "Hnum", "Ldor"]

nymphalid_out_list = ["Mcin", "Dple", "Bany"]

with open(infilename) as infile:
    seq = ''
    for count,line in enumerate(infile):
        if line[0] == ">":
            if count > 1:
                num_chars = len(seq)
                the_as = seq.count("A")
                the_as += seq.count("a")
                the_cs = seq.count("C")
                the_cs += seq.count("c")
                the_gs = seq.count("G")
                the_gs += seq.count("g")
                the_ts = seq.count("T")
                the_ts += seq.count("t")
                the_all = the_as + the_cs + the_gs + the_ts
                if the_all/num_chars >= 0.70:
                    taxon_seq_dict[taxon] = seq
                seq = ''
                taxon = line.strip()[1:]
            else:
                taxon = line.strip()[1:]
        else:
            seq += line.strip()
    num_chars = len(seq)
    the_as = seq.count("A")
    the_as += seq.count("a")
    the_cs = seq.count("C")
    the_cs += seq.count("c")
    the_gs = seq.count("G")
    the_gs += seq.count("g")
    the_ts = seq.count("T")
    the_ts += seq.count("t")
    the_all = the_as + the_cs + the_gs + the_ts
    if the_all/num_chars >= 0.70:
        taxon_seq_dict[taxon] = seq

    for the_taxon in taxon_seq_dict:
        if the_taxon in heliconius_list:
            num_heliconius += 1
        elif "Etal" in the_taxon:
            e_and_a += 1
        elif "Avan" in the_taxon:
            e_and_a += 1
        elif the_taxon in nymphalid_out_list:
            num_nymphalid += 1
    if num_heliconius >= 12 and e_and_a == 2 and num_nymphalid >= 2:
        with open(outfilename, "w") as outfile:
            for the_taxon in taxon_seq_dict:
                #print(len(taxon_seq_dict[the_taxon]))
                outfile.write(">" + the_taxon + "\n")
                outfile.write(taxon_seq_dict[the_taxon] + "\n")
