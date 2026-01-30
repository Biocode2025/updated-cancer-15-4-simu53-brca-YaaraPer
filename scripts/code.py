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



#תוכנית ראשית
BRCA_gene = input("Whether the woman has a familial genetic background and whether she carries a single mutation in BRCA1 or BRCA2, or not: ")
# נהפוך את הקלט לאותיות קטנות כדי להקל על בדיקתו
BRCA_gene = BRCA_gene.lower()





