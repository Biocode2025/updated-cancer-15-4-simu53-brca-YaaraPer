import random
change_genome = ""
seq = "ATGACGGTA"
nucleotide_list = ['T','G','C','A']
if len(seq) == 0:
    print("seq")

rand_num = random.randrange(0,len(seq))
rand_nucleotide = random.choice(nucleotide_list)

nucleotide_list.remove(seq[rand_num])

if seq[rand_num] != rand_nucleotide:
    change_genome = seq[0:rand_num]+ rand_nucleotide + seq[(rand_num+1):]

else:
    nucleotide_list.remove(rand_nucleotide)
    rand_nucleotide = random.choice(nucleotide_list)
    change_genome = seq[0:rand_num]+ rand_nucleotide + seq[(rand_num+1):]