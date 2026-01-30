import random

def Comp_seq(old,new):
  '''
  הפונקציה בודקת כמה הבדלים קיימים בין הרצפים השונים ומחזירה את מספר ההבדלים.
  מקבלת: old,new.
  מחזירה: num_differences.
  '''
  num_differences = 0
 
  for i in range(min(len(old), len(new))):
    if old[i] != new[i]:
      num_differences = num_differences + 1
  return num_differences  
#------------------------------------------------

def Insert_DNA(seq):
  '''
  הפונקציה תכניס במיקום אקראי לרצף ה- DNA של הגן נוקלאוטיד נוסף.
  מקבלת: seq.
  מחזירה: change_genome.
  '''
  nucleotide_list = ['T','G','C','A']
 
  rand_nucleotide = random.choice(nucleotide_list)
  rand_num = random.randrange(0,len(seq))
 
  change_genome = seq[0:rand_num]+ rand_nucleotide + seq[rand_num:]
 
  return change_genome
#------------------------------------------------
 
def Delete_DNA(seq):
  '''
  הפונקציה תחסיר נוקלאוטיד במיקום רנדומאלי.
  מקבלת: seq.
  מחזירה: change_genome.
  '''
  rand_num = random.randrange(0,len(seq))
  rand_nucleotide = seq[rand_num]
 
  change_genome = seq[:rand_num] + seq[rand_num + 1:]
  return change_genome
#------------------------------------------------



#תוכנית ראשית
BRCA_gene = input("Whether the woman has a familial genetic background and whether she carries a single mutation in BRCA1 or BRCA2, or not: ")
# נהפוך את הקלט לאותיות קטנות כדי להקל על בדיקתו
BRCA_gene = BRCA_gene.lower()





