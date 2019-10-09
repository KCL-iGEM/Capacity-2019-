'''
Variable Lists
This will contains call ins for the variables from the csv files and variable lists for funcitons in P2.
'''
import csv
#Data for Disease Table
with open('preset.csv') as csvfileF:
    readcsvF = csv.reader(csvfileF, delimiter=',')
    disease = []
    geneR = []
    gen_sizeR = []
    for rowd in readcsvF:
        dis = rowd[0]
        gen = rowd[1]
        gens = rowd[2]
        disease.append(dis)
        geneR.append(gen)
        gen_sizeR.append(gens)
    disease.pop(0)
    geneR.pop(0)
    gen_sizeR.pop(0)
#Data for radius/ else dimensions
with open('viperdb_clean.csv') as csvfile:
    readcsv = csv.reader(csvfile, delimiter=',')
    inrad = []
    genera = []
    family = []
    T_val = []
    pc_val = []
    genome = []
    total_positions = []
    total_positions2 = []
    for row in readcsv:
        inner_rad = row[6]
        genera_vec = row[2]
        family_vec = row[1]
        T_value = row[4]
        pack_capcity = row[11]
        gen = row[3]
        posi = row[12]
        posi2 = row[12]
        inrad.append(inner_rad)
        genera.append(genera_vec)
        family.append(family_vec)
        T_val.append(T_value)
        genome.append(gen)
        pc_val.append(pack_capcity) #contains all the values for packaging capacity
        total_positions.append(posi)
        total_positions2.append(posi2)
    inrad.pop(0)
    genera.pop(0)
    family.pop(0)
    genome.pop(0)
    T_val.pop(0)
    pc_val.pop(0)
    total_positions.pop(0)
    total_positions2.pop(0)

#protein ratio calculations call ins
 #call in the values for TE, Pex, and Kr - calculation of TE and Comparisons with preex vals for eval selection
with open('tevals.csv') as csvfileT:
    readcsvT = csv.reader(csvfileT, delimiter=',')
    te1 = []
    te2 = []
    te3 = []
    Px1 = []
    Px2 = []
    Px3 = []
    Kr1 = []
    Kr2 = []
    Kr3 = []
    for rowt in readcsvT:
        t1 = rowt[0]
        t2 = rowt[1]
        t3 = rowt[2]
        p1 = rowt[3]
        p2 = rowt[5]
        p3 = rowt[7]
        k1 = rowt[4]
        k2 = rowt[6]
        k3 = rowt[8]
        te1.append(t1)
        te2.append(t2)
        te3.append(t3)
        Px1.append(p1)
        Px2.append(p2)
        Px3.append(p3)
        Kr1.append(k1)
        Kr2.append(k2)
        Kr3.append(k3)
    te1.pop(0)
    te2.pop(0)
    te3.pop(0)
    Px1.pop(0)
    Px2.pop(0)
    Px3.pop(0)
    Kr1.pop(0)
    Kr2.pop(0)
    Kr3.pop(0)
    Kr1 = list(map(float, Kr1))#convert all values to int
    Px1 = list(map(float, Px1))
og_numfacetot = [5, 5, 50]
scal_og = [1, 7.9, 18.6, 28.5, 3.2]
lisdescp = ['P1RBS1', 'P1RBS2', 'P2RBS1', 'P2RBS2', 'P2RBS3']

# create lists for ds versus ss capsid volume evaluations

capsidss = [3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 15, 16, 19, 20, 21, 22, 24, 29, 30, 31, 32, 35, 36, 37, 38, 39, 52, 54,
            57, 59, 60, 61, 62, 63, 64, 65, 67]
used_vectors = [0, 18, 29, 30, 31, 32, 52, 62]  # used in clinical trials - not all are approved for use
used_vectors_once = [0, 18, 29, 52, 62]
parvov_pos = [29, 30, 31, 32]  # mutiples for parvoviridae
approv_vec = [0, 18, 29, 30, 31, 32,
              52]  # this holds all the viral capsids(include repeats) for approved viral vectors
approv_vec_once = [0, 18, 31, 52]  # without repetition of parv vectors
# lists for available capsids

#variables for triangulation number calculations - in global , available to all?
pak_fam = ['Adenoviridae', 'Herpesviridae', 'Parvoviridae', 'Retroviridae']
find_PC_again = [8010000, 32040000, 2403000, 4272000]
sub_ogR = [1500, 960, 60, 60]
T_famR = [25, 16, 1, 1]
scalR = [91.16052458, 53.359577, 17.652715, 2.541513155]
vol_outerradR = [478829133.8, 1121104558, 27816758.06, 7119751.55]


#creating a summary table - display all T val, and potential readjustrments
#final summary table - print for all conditions?
pak_famB = ['Adenoviridae', 'Herpesviridae', 'Parvoviridae', 'Retroviridae']
sub_ogRB = [1500, 960, 60, 60]
T_famRB = [25, 16, 1, 1]
scalRB = [91.16052458, 53.359577, 17.652715, 2.541513155]
vol_outerradRB = [478829133.8, 1121104558, 27816758.06, 7119751.55]

#CAPSID stability
stability_vir = [-0.42592213, -1.013525132, -0.027586207]
energies_vir = [-415.7, -383.1125, -3.2]
outer_diam = [976, 378, 116]
pk_lis = ['Adenoviridae', 'Parvoviridae', 'Retroviridae']






