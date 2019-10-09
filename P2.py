'''
Functions for CapsidOptimizer and CapsidBuilder Calculations
'''
from P1 import*
#create table for disease list
def table_presets():
    print("IS THERE A VIRAL VECTOR THAT CAN DELIVERY YOUR GENE?\n\n")
    # create summary table for the presets
    print("\t" * 18 + "SUMMARY RARE GENETIC DISEASES")
    titles = ['Disease', 'Gene Name', 'Gene Length(nt)']
    dataF = [titles] + list(zip(disease, geneR, gen_sizeR))
    for i, d in enumerate(dataF):
        line = ''.join(str(x).ljust(58) for x in d)
        print("\t" + line)
        if i == 0:
            print('\t' + '-' * len(line) + '\n')
def ask_user():
    global gene
    gene = input("\nPlease enter a gene length from the table above(nt): ")
    return gene

def gene_vol():
    global gen_vol, gen_volss
    gen_vol = int(gene)
    gen_vol = (1068 * gen_vol) + 200  # this convert to a ds clindrical volume
    gen_volss = (534 *int(gene)) + 200
    return gen_vol, gen_volss

def calc_sphericalVol():
    import math
    global inrad
    global capsid_final
    inrad = list(map(int, inrad))  # convert all the inner capsid radi to integer (read in as string)
    # calcualte the spherical volume of the capsid shells (in viperdb) - from inner radius
    list_invol = []
    for inrad_vol in inrad:
        inrad_vol = 4.0 / 3.0 * math.pi * inrad_vol ** 3
        list_invol.append(inrad_vol)
    capsid_final = list_invol
    return capsid_final


def list_capsidVOL_ds():
    global capss_vol, new_capsid_final
    capss_vol = []
    for sids in capsidss:
        lid = capsid_final[sids]  # will find the volume at this position in the total list
        capss_vol.append(lid)  # contains all the volumes of the capsid with ss genome type
    new_capsid_final = capsid_final.copy()
    for lose in capss_vol:
        if lose in new_capsid_final:
            new_capsid_final.remove(lose)  # this is the list of only ds capsid volumes
    return capss_vol, new_capsid_final

def compare_genomeandcapsidVol():
    global counter, dec1, dec2, c1_lisgen, c1_lisvol, c1_posav, c2_lisgen, c2_lisvol
    c1_lisgen = []
    c1_lisvol = []
    c1_posav = []
    c2_lisgen = []
    c2_lisvol = []
    dec1 = False
    dec2 = False
    counter = 0
    for new_cap in range(68):
        var_vol = capsid_final[new_cap]
        if var_vol in new_capsid_final:  # is it a capsid with ds genome?
            if gen_vol < var_vol:
                counter += 1
                c1_lisgen.append(genera[new_cap])  # add the name of available genera
                c1_lisvol.append(var_vol)  # add the volume of the avail capsid
                c1_posav.append(new_cap)  # this will hold the all the available capsids
                dec1 = True
                yield dec1, 'dec1', counter, c1_posav
            elif gen_vol > var_vol:  # capsid can not hold the gene based on volume alone
                c2_lisgen.append(genera[new_cap])  # hold the genera that do not fit the gene
                c2_lisvol.append(var_vol)
                dec2 = True
                yield dec2, 'dec2', counter, 'dec0', dec1
        elif var_vol in capss_vol:
            if gen_volss < var_vol:  # capsid can hold the gene - based on volume alone
                counter += 1
                c1_lisgen.append(genera[new_cap])  # add the name of available genera
                c1_lisvol.append(var_vol)  # add the volume of the avail capsid
                c1_posav.append(new_cap)  # this will hold the all the available capsids
                dec1 = True
                yield dec1, 'dec1', counter, c1_posav
            elif gen_volss > new_cap:  # capsid can not hold the gene based on volume alone
                c2_lisgen.append(genera[new_cap])  # hold the genera that do not fit the gene
                c2_lisvol.append(var_vol)
                dec2 = True
                yield dec2, 'dec2', counter, 'dec0', dec1

def nocapsidsVol():
    global conv_sphericalV
    pi = 3.14159265359
    print("\nYour gene can not fit inside any available viral capsids!\n *determined from inner capsid volume only")
    conv_sphericalV = (2 / 3) * gen_vol  # convert the cylindrical volume to a spherical volume
    est_inrad = ((3 * conv_sphericalV) / (4 * pi)) ** (1. / 3.)
    print("\n\nThe minimum packaging requirements for your gene are as listed:")
    print("---------------Theoretical Viral Capsid Dimensions--------------")
    print(" Minimum Capsid  (A): " + str(int(conv_sphericalV)))
    print(" New inner radius: " + str(int(est_inrad)) + " A")
    print("------------------------------------------------------------------")
    print("\nThe capsid(s) unavailable for your gene include: ")
    titles = ['Genera', 'Volume (A)']
    data1 = [titles] + list(zip(c2_lisgen, c2_lisvol))
    for i, d in enumerate(data1):
        line = ''.join(str(x).ljust(30) for x in d)
        print("\t" + line)
        if i == 0:
            print('\t' + '-' * len(line) + '\n')
def allcapsidsVOL():
    global conv_sphericalV
    pi = 3.14159265359
    conv_sphericalV = (2 / 3) * gen_vol  # convert the cylindrical volume to a spherical volume
    est_inrad = ((3 * conv_sphericalV) / (4 * pi)) ** (1. / 3.)
    print("\n\nThe minimum packaging requirements for your gene are as listed:")
    print("---------------Theoretical Viral Capsid Dimensions--------------")
    print(" Minimum Capsid  (A): " + str(int(conv_sphericalV)))
    print(" New inner radius: " + str(int(est_inrad)) + " A")
    print("------------------------------------------------------------------")
    print("\nYour gene fits in " + str(counter) + " out of 68 available viral capsids based on capsid volume alone..")
    print("\nThe viral capsid(s) available for your gene are as listed:")
    titles = ['Genera', 'Volume (A)']
    data2 = [titles] + list(zip(c1_lisgen, c1_lisvol))
    for i, d in enumerate(data2):
        line = ''.join(str(x).ljust(30) for x in d)
        print("\t" + line)
        if i == 0:
            print('\t' + '-' * len(line) + '\n')
def somecapsidsVOL():
    global conv_sphericalV
    pi = 3.14159265359
    print("\n\nThe minimum packaging requirements for your gene are as listed:")
    conv_sphericalV = (2 / 3) * gen_vol  # convert the cylindrical volume to a spherical volume
    est_inrad = ((3 * conv_sphericalV) / (4 * pi)) ** (1. / 3.)
    print("\n\nThe minimum packaging requirements for your gene are as listed:")
    print("---------------Theoretical Viral Capsid Dimensions--------------")
    print(" Minimum Capsid  (A): " + str(int(conv_sphericalV)))
    print(" New inner radius: " + str(int(est_inrad)) + " A")
    print("------------------------------------------------------------------")
    print("Your gene only fits in " + str(counter) + " out of 68 available viral capsids based on capsid volume alone.")
    print("\nThe viral capsid(s) available for your gene are as listed:")
    titles = ['Genera', 'Volume (A)']
    data2 = [titles] + list(zip(c1_lisgen, c1_lisvol))
    for i, d in enumerate(data2):
        line = ''.join(str(x).ljust(30) for x in d)
        print("\t" + line)
        if i == 0:
            print('\t' + '-' * len(line) + '\n')
    print("\nThe remaining capsid(s) unavailable for use include: ")
    titles = ['Genera', 'Volume (A)']
    data1 = [titles] + list(zip(c2_lisgen, c2_lisvol))
    for i, d in enumerate(data1):
        line = ''.join(str(x).ljust(30) for x in d)
        print("\t" + line)
        if i == 0:
            print('\t' + '-' * len(line) + '\n')


def approvlisting():
    global fam_lisF, fam_lisF2, gen_listT, fam_listT, capvol_ap, capvol_apR, famvol_apR
    # FINDING APPROVED VIRAL CAPSID
    fam_lisF = []  # holds viral vectors that are being used
    fam_lisF2 = []  # holds vectors that have approved therapies

    gen_listT = []
    fam_listT = []
    capvol_ap = []  # hold the volumes of the used vectors (not approved for therapeutics)

    # readj table variables
    capvol_apR = []  # holds the volumes of only the approved viral vectors
    famvol_apR = []  # holds the family name of the approved viral vectors

    for vol2 in approv_vec:  # release values for positions of approved vectors
        volad2 = capsid_final[vol2]  # will find the volumes from the positions
        capvol_apR.append(volad2)  # this will hold all the capsid volume of approved vectors

    for famVec in approv_vec:
        famVar = family[famVec]
        famvol_apR.append(famVar)  # hold all the families of the approved viral vectors

    for vol in used_vectors:
        volad = capsid_final[vol]
        capvol_ap.append(volad)

    for item in used_vectors_once:
        fam_lisI = family[item]
        fam_lisF.append(fam_lisI)

    for item3 in used_vectors:
        gen_lis_T = genera[item3]
        gen_listT.append(gen_lis_T)  # hold genera for viral families that are being used

    for item4 in used_vectors:
        fam_lis_T = family[item4]
        fam_listT.append(fam_lis_T)  # hold family of viruses being used

    for item2 in approv_vec_once:
        fam_lisI2 = family[item2]
        fam_lisF2.append(fam_lisI2)

    return fam_lisF, fam_lisF2, gen_listT, fam_listT, capvol_ap

def findapprovedcapsids():
    global gen_approvlis, fam_approvlis, pacap_approvlis, record_pos
    gen_approvlis = []  # THIS WILL HOLD ONLY ONES TESTED AGAINST THE GENE INPUT
    fam_approvlis = []
    pacap_approvlis = []
    record_pos = []
    cond_2 = False
    for pos in c1_posav:
        if pos in approv_vec_once:
            cond_2 = True
            genap = genera[pos]
            famapp = family[pos]
            pacapp = pc_val[pos]
            record_pos.append(pos)
            gen_approvlis.append(genap)
            fam_approvlis.append(famapp)
            pacap_approvlis.append(pacapp)
    return cond_2
def findaprovTABLE():
    print("\n\t\t\t\t\tSUMMARY TABLE OF APPROVED VIRAL VECTORS WITH LARGER CAPSID VOLUME THAN YOUR GENE\n")
    titles = ['Family', 'Genera', 'Volume (A)']
    data11 = [titles] + list(zip(fam_approvlis, gen_approvlis, capvol_apR))
    for i, d in enumerate(data11):
        line = ''.join(str(x).ljust(26) for x in d)
        print("\t" + line)
        if i == 0:
            print('\t' + '-' * len(line) + '\n')

def someavailableSUM():
    find_sml_vol = []
    for recP in record_pos:
        find_vol = capsid_final[recP]
        find_sml_vol.append(find_vol)
    smallest_vol = min(find_sml_vol)
    dest = find_sml_vol.index(smallest_vol)
    dest2 = record_pos[dest]
    print("\n*Of these approved viral vectors, the one with the smallest CAPSID VOLUME is the " + family[
        dest2] + " viral vector. Which has a volume of " + str(int(smallest_vol)) + " A")
    print("-" * 180)

    # find the positions of the capsids that are unavailable
    remainder = []
    for fin in c1_posav:
        fin = str(fin)
        if fin in total_positions:  # this will identify the string duplicate in both lists
            fin = int(fin)
            remainderval = total_positions[fin]  # collect these positions in a variable
            remainder.append(remainderval)
    for bin in remainder:
        if bin in total_positions:
            total_positions.remove(
                bin)  # total position list now changed to only include the positions of unaval capsids
    fam_unavl = []
    gen_unavl = []
    pac_unavl = []
    vol_unavl = []
    for sin in approv_vec_once:
        sin = str(sin)
        if sin in total_positions:
            sin = int(sin)
            nameF_unaval = family[sin]
            nameG_unaval = genera[sin]
            nameP_unaval = pc_val[sin]
            nameV_unaval = capsid_final[sin]
            fam_unavl.append(
                nameF_unaval)  # these are the lists of information for all the remaining unavailble capsids
            gen_unavl.append(nameG_unaval)
            pac_unavl.append(nameP_unaval)
            vol_unavl.append(nameV_unaval)
    # create a list of readjustments
    dif1_volume = []  # contains the readjustments for the approved viral capsids
    ad1_volume = []
    for num1 in vol_unavl:
        if num1 in new_capsid_final:
            dif1 = float(gen_vol) - float(num1)
            dif1_volume.append(dif1)
            for num4, num5 in zip(dif1_volume,
                                  vol_unavl):  # unpack, both difference list and the approved vectors - add final will all be the same (DO NEED THIS LIST?)
                add_on1 = num4 + num5  # this will add the difference to the list of approved capsid volumes
                ad1_volume.append(add_on1)  # contain the readjusted volumes in a new list
        if num1 in capss_vol:  # this find if any of the unavl capsids are inside the list of ss capsid -
            dif2 = float(gen_volss) - float(num1)
            dif1_volume.append(dif2)
            for num6, num7 in zip(dif1_volume, vol_unavl):
                add_on2 = num6 + num7
                ad1_volume.append(add_on2)

        # table for readjusted values --- THIS IS BASED ONLY ON CAPSID VOLUME
    print("\n\n\t\t\t\t\t\t\t_READJUSTMENTS OF APPROVED VIRAL *CAPSID VOLUMES_\n")
    titles = ['Family', 'Genera', 'Capsid Volume', 'Required Adjustment', 'New Capsid Volume']
    data = [titles] + list(zip(fam_unavl, gen_unavl, vol_unavl, dif1_volume, ad1_volume))
    for i, d in enumerate(data):
        line = ''.join(str(x).ljust(25) for x in d)
        print('\t\t\t\t' + line)
        if i == 0:
            print('\t\t\t' + '-' * len(line))
    print("\nThe minimum capsid volume requirements for your gene is: " + str(int(conv_sphericalV)) + " A")
    print("The genome volume *(ds measurements) for your gene is: " + str(gen_vol) + "A")
    print("The genome volume *(ss measurements) for your gene is " + str(gen_volss) + "A")

def noavailableSUM():
    dif_volume = []  # contains the readjustments for the approved viral capsids
    ad_volume = []
    print('\nYOUR GENE CAN NOT BE DELIVERED BY ANY APPROVED VIRAL VECTORS.'
          '\nTHEIR INNER CAPSID VOLUME MUST BE INCREASED TO ACCOMODATE YOUR GENE')
    for num in capvol_apR:  # iterate over the values of the approved vectors - no selection needed, assumption all unavailable
        dif = float(gen_vol) - float(num)  # how much need to increase? -
        dif_volume.append(
            dif)  # holds the readjusted capsids values - how must need to increase the capsids based on the genome vol and approved vectors

    for num2, num3 in zip(dif_volume,
                          capvol_apR):  # unpack, both difference list and the approved vectors - add final will all be the same (DO NEED THIS LIST?)
        add_on = num2 + num3  # this will add the difference to the list of approved capsid volumes
        ad_volume.append(add_on)  # contain the readjusted volumes in a new list

    # table for readjusted values --- THIS IS BASED ONLY ON CAPSID VOLUME
    print("\n\n\t\t\t\t\t\t\t_READJUSTMENTS OF REMAINING APPROVED VIRAL *CAPSID VOLUMES_\n")
    titles = ['Family', 'Capsid Volume', 'Required Adjustment', 'New Capsid Volume']
    data = [titles] + list(zip(famvol_apR, capvol_apR, dif_volume, ad_volume))
    for i, d in enumerate(data):
        line = ''.join(str(x).ljust(25) for x in d)
        print('\t\t\t\t' + line)
        if i == 0:
            print('\t\t\t' + '-' * len(line))
    print("\nThe minimum capsid volume requirements for your gene is: " + str(conv_sphericalV) + " A")
    print("The genome volume *(ds measurements) for your gene is: " + str(gen_vol) + "A")

def convertpctoVOL():
    global capsVlis
    pak_capV = [7500, 30000, 4500, 8000]  # no repeats for parvovirdae
    capsVlis = []
    for caps in pak_capV:
        capsV = 1068 * caps
        capsVlis.append(capsV)  # THIS WILL HOLD ALL THE packing volumes
    target = [2403000, 4272000]
    pos_change = [2, 3]
    for x, y in zip(pos_change, target):
        capsVlis[x] = y
    return capsVlis

def sortandfindPCvals():
    global pp1, pp2, pp1_V, pp2_V
    count_pack = -1  # this will record the names of positions in both conditions
    pp1 = []  # will hold family names of the available capsids
    pp2 = []  # will hold the family names of the unavailable capsids ---- both are based on packing volumes(cylindrical)
    pp1_V = []
    pp2_V = []  # this will hold the volumes of the unavailable capsids
    pak_fam = ['Adenoviridae', 'Herpesviridae', 'Parvoviridae', 'Retroviridae']
    T_fam = [25, 16, 1, 1]  # holds the T values of the approved capsids
    counter_packingV = 0  # this used for output statements later
    case1 = False
    case2 = False
    re1 = 0  # this will record the positions of the ones that are approved
    re1store = []  # this will store all the available positions
    # Unpack volumes of the calculated volumes within a loop
    for find1 in capsVlis[0:2]:  # between only the ds capsid volumes
        if int(gen_vol) < int(find1):  # THIS CAN BE USED FOR THE GENE - APPROVED PACKING CAPACITY - this compares a ds gene to capsV_lis - should change according to same genome type
            counter_packingV += 1  # this is used later for the output list
            count_pack += 1  # record the position inside the list ---use this to find the family namesw
            new_pack2 = pak_fam[count_pack]  # - use the position to find the name of the capsid
            pp1.append(new_pack2)  # collect all the names of the capsids larger than cargo ---can be used
            pp1_V.append(find1)  # this will collect only available volumes (calc of pack)
            case1 = True
            re1 += 1
            re1store.append(
                re1)  # hold onto all the positions that are added -- this order will be printed in the list below
        if int(gen_vol) > int(find1):  # THIS CAN NOT BE USED FOR THE GENE ---APPROVED PACKING CAPACITY
            count_pack += 1
            new_pack = pak_fam[count_pack]
            pp2.append(new_pack)  # collect all the capsids smaller than cargo -- can not be used
            pp2_V.append(find1)
            case2 = True
        if count_pack == 1:  # when the loop is finished evaluating all the ds genome type then go to the next two ss genome types
            for find2 in capsVlis[2:5]:
                if int(gen_volss) < find2:
                    counter_packingV += 1
                    count_pack += 1
                    new_pack3 = pak_fam[count_pack]
                    pp1.append(new_pack3)  # collect the names of the ss capsids
                    pp1_V.append(find2)  # collect the volumes of ss capsid volumes
                    case1 = True
                    re1 += 1
                    re1store.append(re1)
                if int(gen_volss) > find2:
                    count_pack += 1
                    new_pack4 = pak_fam[count_pack]
                    pp2.append(new_pack4)  # add the capsid name to the pp2 list
                    pp2_V.append(find2)  # this will collect the capsids volumes that are unavailable
                    case2 = True
    if counter_packingV == 0:
        print("These are the approved viral vectors with INSUFFICIENT PACKAGING CAPACITIES THAT can not DELIVER YOUR GENE: " + ", ".join(pp2))
    elif case1 == True and case2 == False:
        print("\nThese are approved viral vectors with SUFFICIENT PACKAGING CAPACITIES THAT can DELIVER YOUR GENE: " + ", ".join(pp1))
    elif case1 and case2 == True:
        print("\nThese are the approved viral vectors with SUFFICIENT PACKAGING CAPACITIES THAT can DELIVER YOUR GENE: " + ", ".join(pp1))
        print("These are the approved viral vectors with INSUFFICIENT PACKAGING CAPACITIES THAT can not DELIVER YOUR GENE: " + ", ".join(pp2))

def addposofFOUNDPC():
    global add_pos1, add_pos2
    add_pos1 = []
    if 'Adenoviridae' in pp1:
        add_pos1.append(7500)
    if 'Herpesviridae' in pp1:
        add_pos1.append(30000)
    if 'Parvoviridae' in pp1:
        add_pos1.append(4500)
    if 'Retroviridae' in pp1:
        add_pos1.append(8000)
    add_pos2 = []
    if 'Adenoviridae' in pp2:
        add_pos2.append(7500)
    if 'Herpesviridae' in pp2:
        add_pos2.append(30000)
    if 'Parvoviridae' in pp2:
        add_pos2.append(4500)
    if 'Retroviridae' in pp2:
        add_pos2.append(8000)
def havePCatleastone():
    global name_cap, smallest_volPC
    # CREATE A SUMMARY TABLE OF ONLY THE CAPSID WITH SUFFICIENT PACKAGING CAPACITY
    print("\n\tSummary Table of Viral Vectors that CAN Package your Gene\n")
    titles = ['Family', 'Packaging Capacity(bp)', 'Packaging Volume (A)']
    data11 = [titles] + list(zip(pp1, add_pos1, pp1_V))
    for i, d in enumerate(data11):
        line = ''.join(str(x).ljust(26) for x in d)
        print("\t" + line)
        if i == 0:
            print('\t' + '-' * len(line) + '\n')
    # within the pack table, find the smallest capsid ---- there must be at least one entry inside this table - this will be the selected capsid
    smallest_volPC = min(pp1_V)  # this will hold the smallest volume
    name_cap_pos = pp1_V.index(smallest_volPC)
    name_cap = pp1[name_cap_pos]
    print("\n*Of these approved viral vectors, the one with the smallest packaging capacity is " + str(
        name_cap) + " viral vector. Which has a volume of " + str(int(smallest_volPC)) + " A")
    print("\n" + "-" * 180)
    # find the genome type of this selected capsid
    if smallest_volPC == 8010000:
        gentype3 = gen_vol
    if smallest_volPC == 32040000:
        gentype3 = gen_vol
    if smallest_volPC == 2403000:
        gentype3 = gen_volss
    if smallest_volPC == 4272000:
        gentype3 = gen_volss

    # CREATE A SUMMARY TABLE OF READJUSTMENTS FOR THIS CONDITION ALSO -- USE PP2 VARIABLE
    adj_pcVf2 = []  # hold the final values (all the same) of the minimum cylindrical volume of the selected gene
    adj_pcVi2 = []  # hold the difference - required readjustment
    # must evaluate the seperate types of genomes individually - find the ds genome sizes within the pp2 list first
    for java in pp2_V:  # begin to iterate over the pp2 list, then create seperate evaluation conditions
        if int(java) == 8010000:
            diff_pcVee2 = gen_vol - int(java)  # only ds nucleic acid
            add_pcVee2 = int(java) + diff_pcVee2
            adj_pcVi2.append(diff_pcVee2)
            adj_pcVf2.append(add_pcVee2)
        if int(java) == 32040000:
            diff_pcVee4 = gen_vol - int(java)  # only ds nucleic acid
            add_pcVee4 = int(java) + diff_pcVee4
            adj_pcVi2.append(diff_pcVee4)
            adj_pcVf2.append(add_pcVee4)
        if int(java) == 2403000:
            diff_pcVee3 = gen_volss - int(java)
            add_pcVee3 = int(java) + diff_pcVee3
            adj_pcVi2.append(diff_pcVee3)
            adj_pcVf2.append(add_pcVee3)
        if int(java) == 4272000:
            diff_pcVee5 = gen_volss - int(java)
            add_pcVee5 = int(java) + diff_pcVee5
            adj_pcVi2.append(diff_pcVee5)
            adj_pcVf2.append(add_pcVee5)
    # reajustment table for when there are some* capsids with a sufficient packaging capcity
    print("\n\n\t\t\t\t\t\t\t_READJUSTMENTS OF APPROVED VIRAL PACKAGING CAPACITY_\n")
    titles = ['Family', 'Packaging Capacity(A)', 'Required Adjustment(A)', 'New Packaging Capacity(A)']
    data = [titles] + list(zip(pp2, pp2_V, adj_pcVi2, adj_pcVf2))
    for i, d in enumerate(data):
        line = ''.join(str(x).ljust(25) for x in d)
        print('\t\t\t\t' + line)
        if i == 0:
            print('\t\t\t' + '-' * len(line))

    low_red2 = min(adj_pcVi2)
    name_capL3 = adj_pcVi2.index(low_red2)  # holds position in a list of 4
    name_capL4 = pp2[name_capL3]
    print("\n*Of these approved vectors, the one with the smallest packaging capacity readjustment is the " + str(
        name_capL4) + " viral vector. Where it must increase it's packaging volume by: " + str(
        int(low_red2)) + " A")
    # create a statement for determining the minimum packaging capcity based on the genome type
    dif1 = int(gen_vol) - int(low_red2)  # variable for either adeno or herp
    dif2 = int(gen_volss) - int(low_red2)  # variable for either parv or retro
    if dif1 == 32040000:
        gentype = gen_vol  # for the selected capsid, it will hold ds nucleic acid
        print("\nThe minimum capsid volume requirements for your gene is: " + str(gentype) + " A.")
        print("\n" + "-" * 180)
    if dif1 == 8010000:
        gentype = gen_vol  # for the selected capsid, it will hold ds nucleic acid
        print("\nThe minimum capsid volume requirements for your gene is: " + str(gentype) + " A.")
        print("\n" + "-" * 180)
    if dif2 == 2403000:
        gentype = gen_volss  # for the selected capsid, it will hold ss nucleic acid genome
        print("\nThe minimum capsid volume requirements for your gene is: " + str(gentype) + " A.")
        print("\n" + "-" * 180)
    if dif2 == 4272000:
        gentype = gen_volss  # for the selected capsid, it will hold ss nucleic acid genome
        print("\nThe minimum capsid volume requirements for your gene is: " + str(gentype) + " A.")
        print("\n" + "-" * 180)
    low_redVol2 = pp2_V[name_capL3]  # this holds the volume of the selected capsid with the smallest readjustment

    # this statement will always be printed for this section - bc there will always be at least one approved capsid in the list
    print(
        "\t\t\t| THE MOST APPROPRIATE APPROVED VIRAL CAPSID WITH AN SUFFICIENT PACKAGING CAPACITY RECOMMENDED FOR THE DELIVERY OF YOUR GENE IS THE: " + str(
            name_cap) + " viral vector |")  # should be the smallest capsid within the list

def noPCone():
    global name_capL1, low_redVol
    # CREATE A SUMMARY TABLE OF ONLY THE CAPSID WITH SUFFICIENT PACKAGING CAPACITY
    print("\n\tSummary Table of Viral Vectors that CAN Package your Gene\n")
    titles = ['Family', 'Packaging Capacity(bp)', 'Usable Volume (A)']
    data11 = [titles] + list(zip(pp1, add_pos1, pp1_V))
    for i, d in enumerate(data11):
        line = ''.join(str(x).ljust(26) for x in d)
        print("\t" + line)
        if i == 0:
            print('\t' + '-' * len(line) + '\n')
    print("THERE ARE NO APPROVED VIRAL VECTORS THAT HAVE A LARGE ENOUGH PACKAGING CAPACITY TO DELIVERY YOUR GENE!")
    # distinguish between the genome type of the capsid, whether it is ds or ss type genome to find the appropriate readjustment values
    adj_pcVf = []
    adj_pcVi = []
    for vee in capsVlis:
        if int(vee) == 8010000:
            diff_pcVee6 = gen_vol - int(vee)  # only ds nucleic acid
            add_pcVee6 = int(vee) + diff_pcVee6
            adj_pcVi.append(diff_pcVee6)
            adj_pcVf.append(add_pcVee6)
        if int(vee) == 32040000:
            diff_pcVee7 = gen_vol - int(vee)  # only ds nucleic acid
            add_pcVee7 = int(vee) + diff_pcVee7
            adj_pcVi.append(diff_pcVee7)
            adj_pcVf.append(add_pcVee7)
        if int(vee) == 2403000:
            diff_pcVee8 = gen_volss - int(vee)
            add_pcVee8 = int(vee) + diff_pcVee8
            adj_pcVi.append(diff_pcVee8)
            adj_pcVf.append(add_pcVee8)
        if int(vee) == 4272000:
            diff_pcVee9 = gen_volss - int(vee)
            add_pcVee9 = int(vee) + diff_pcVee9
            adj_pcVi.append(diff_pcVee9)
            adj_pcVf.append(add_pcVee9)

    print("\n\n\t\t\t\t\t\t\t_READJUSTMENTS OF APPROVED VIRAL PACKAGING CAPACITY_\n")
    titles = ['Family', 'Packaging Capacity(A)', 'Required Adjustment(A)', 'New Packaging Capacity(A)']
    data = [titles] + list(zip(pp2, capsVlis, adj_pcVi, adj_pcVf))
    for i, d in enumerate(data):
        line = ''.join(str(x).ljust(25) for x in d)
        print('\t\t\t\t' + line)
        if i == 0:
            print('\t\t\t' + '-' * len(line))
    # find the one requires the smallest readjustment
    low_red = min(adj_pcVi)
    name_capL2 = adj_pcVi.index(low_red)  # holds position in a list of 4
    name_capL1 = pp2[name_capL2]
    print("\n*Of these approved viral vectors, the one with the smallest packaging capacity readjustment is " + str(
        name_capL1) + " viral vector. Where it must increase it's packaging volume by: " + str(int(low_red)) + " A")
    # create a statement for determining the minimum packaging capcity based on the genome type
    dif3 = int(gen_vol) - int(low_red)  # variable for either adeno or herp
    dif4 = int(gen_volss) - int(low_red)  # variable for either parv or retro
    if dif3 == 32040000:
        gentype2 = gen_vol  # for the selected capsid, it will hold ds nucleic acid
        print("\nThe minimum capsid volume requirements for your gene is: " + str(gentype2) + " A.")
        print("\n" + "-" * 180)
    if dif3 == 8010000:
        gentype2 = gen_vol  # for the selected capsid, it will hold ds nucleic acid
        print("\nThe minimum capsid volume requirements for your gene is: " + str(gentype2) + " A.")
        print("\n" + "-" * 180)
    if dif4 == 2403000:
        gentype2 = gen_volss  # for the selected capsid, it will hold ss nucleic acid genome
        print("\nThe minimum capsid volume requirements for your gene is: " + str(gentype2) + " A.")
        print("\n" + "-" * 180)
    if dif4 == 4272000:
        gentype2 = gen_volss  # for the selected capsid, it will hold ss nucleic acid genome
        print("\nThe minimum capsid volume requirements for your gene is: " + str(gentype2) + " A.")
        print("\n" + "-" * 180)
    print(
        "\t" * 6 + "| THE MOST APPROPRIATE APPROVED VIRAL CAPSID RECOMMENDED FOR THE DELIVERY OF YOUR GENE IS THE: " + str(
            name_capL1) + " viral vector |")
    # then find its original volume for later!!!
    low_redVol = capsVlis[name_capL2]  # hols the volume of the selected capsid above -- use this for t value!

def availablePCTVAL():
    vol_outerdes = []
    find_ogvol = []
    store_faceR = []
    gen_typer = []
    global decide_rbs1, trip, original_voly, rmR, volume_desiredR, volume_ic1R, oldT, oldprot, rmR, rm_new1R
    decide_rbs1 = True
    find_po = pak_fam.index(str(name_cap))  # find the position of the selected capsid in the name list
    trip = find_po  # equivalent matching pos in other lists items
    if trip < 2:
        original_voly = vol_outerradR[trip]
        find_ogvol.append(original_voly)
        genometypeR = gen_vol
        gen_typer.append(genometypeR)
        gen_spR = genometypeR * (2 / 3)  # convert to the volume of a sphere
        gen_icosR = gen_spR * (0.983631649)  # convert to a icosahedron
        # convert this to a capsid with a new overall outer radius(this will be different for each capsid - scaling value)
        gen_changR = gen_icosR * scalR[
            trip]  # for the original conversion between inner and outer capsid volumes for each virus
        volume_desiredR = gen_changR  # the volume of our gene will determines the minimum required/target volume of our capsid
        vol_outerdes.append(volume_desiredR)
        pi = 3.14159265359
        # find the area of the triangular face
        rmR = 0.8090169944 * ((volume_desiredR / 2.181694991) ** (1. / 3.))  # find a radius for the target capsid
        side_lengthR = rmR / (0.80901699)
        total_areaR = (8.6602540378) * ((side_lengthR) ** 2)
        area_faceR = total_areaR / 20  # single face
        store_faceR.append(area_faceR)
        print("The area of a single triangular face from your required icosahedral capsid is: " + str(
            int(area_faceR)) + " A^2.")
        print(" Which had a volume of: " + str(gen_vol) + " A^3")
    elif trip >= 2:
        original_voly = vol_outerradR[trip]
        find_ogvol.append(original_voly)
        genometypeR1 = gen_volss
        gen_typer.append(genometypeR1)
        gen_spR = genometypeR1 * (2 / 3)  # convert to the volume of a sphere
        gen_icosR = gen_spR * (0.983631649)  # convert to a icosahedron
        # convert this to a capsid with a new overall outer radius(this will be different for each capsid - scaling value)
        gen_changR = gen_icosR * scalR[
            trip]  # for the original conversion between inner and outer capsid volumes for each virus
        volume_desiredR = gen_changR  # the volume of our gene will determines the minimum required/target volume of our capsid
        vol_outerdes.append(volume_desiredR)
        pi = 3.14159265359
        # find the area of the triangular face
        rmR = 0.8090169944 * ((volume_desiredR / 2.181694991) ** (1. / 3.))  # find a radius for the target capsid
        side_lengthR = rmR / (0.80901699)
        total_areaR = (8.6602540378) * ((side_lengthR) ** 2)
        area_faceR = total_areaR / 20  # single face
        store_faceR.append(area_faceR)
        print("The area of a single triangular face from your required icosahedral capsid is: " + str(
            int(area_faceR)) + " A^2.")
        print("Which had a volume of: " + str(gen_volss) + " A^3")

    # find the volume of the original capsid
    volume_ic1 = smallest_volPC
    seeker = find_PC_again.index(volume_ic1)
    volume_ic1R = vol_outerradR[int(seeker)]
    rm_new1R = 0.8090169944 * ((volume_ic1R / 2.181694991) ** (1. / 3.))
    side_length1R = rm_new1R / (0.80901699)
    total_area1R = (8.6602540378) * ((side_length1R) ** 2)
    area_face1R = total_area1R / 20  # area of the original side face
    print("The area of a single triangular face of the " + str(name_cap) + "viral vector is: "
          + str(int(area_face1R)) + " A^2." + " Which had a volume of: " + str(int(smallest_volPC)) + " A")
    # find the new T value
    store_newTR = []
    store_newprot = []
    store_enlg = []
    oldT = round(int(T_famR[trip]))
    oldprot = round((int(sub_ogR[trip])))
    for x in store_faceR:
        enlargeR = int(x) / area_face1R
        store_enlg.append(enlargeR)
        newTR = round(int(oldT) * enlargeR)
        store_newTR.append(newTR)  # hold all the new t values
        num_protein1R = 60 * (newTR)
        store_newprot.append(num_protein1R)
        # this will store the new number of protein required
    enlargRatio = store_enlg[0]
    new_amountofprotein = store_newprot[0]
    more_protein = int(new_amountofprotein) - int(oldprot)
    pick_vol = vol_outerdes[0]
    for ttr in store_newTR:
        if enlargRatio < 1:
            print("\nThe enlargement ratio between your required capsid and selected capsid is: " + str(enlargRatio) + "\nThis means the new capsid has the same triangulation number of the original capsid which was " + str(
                oldT))
            print( "\n*From this result, we can confirm that this viral capsid is most appropriate choice for the delivery of your gene and does not require further readjustments.")
            print("\n\n\t\t\t\t\t\t\t_Summary of Viral Capsid Optimization_\n")
            titles = ['Selected Capsid', 'Capsid Volume(A)', 'T Value(A)',
                      'No. Subunits']
            data = [name_cap, volume_ic1R, oldT, oldprot]
            print("\t" + "|  " + "    |   ".join(titles) + " |")
            print("\t" + "-" * 87)
            print("\t   " + "\t\t\t\t".join(map(str, data)))
        elif enlargRatio >= 1:
            print("\nThe enlargement ratio between your required capsid and selected capsid is: " + str(
                enlargRatio) + "\nThis means the new capsid has a triangulation number of " + str(ttr))
            print("To deliver your gene, the original " + str(
                name_cap) + "viral vector must be increased in size by a factor of " + str(enlargRatio))
            print("In order to acheive this, the number of proteins that make up this viral capsid will need to increase by " + str(
                    round(more_protein)))
            print("\n\n\t\t\t\t\t\t\t_Summary of Viral Capsid Optimization_\n")
            titles = ['Selected Capsid', ' Original Capsid Volume(A)', 'Original T Value(A)',
                      'Original No. Subunits', 'New Capsid Volume(A)', 'New T Value(A)', 'New No. Subunits']
            data = [name_cap, volume_ic1R, oldT, oldprot, pick_vol, ttr, new_amountofprotein]
            print("|" + "  | ".join(titles) + "|")
            print("-" * 180)
            print("\t\t\t\t\t".join(map(str, data)))
    return decide_rbs1
def NOavailablePCTVAL():
    vol_outerdes = []
    find_ogvol = []
    store_faceR = []
    gen_typer = []
    global decide_rbs2, trip, original_voly, rmR, volume_desiredR, volume_ic1R, oldT, oldprot, ttr, new_amountofprotein, gen_vol, rmR, rm_new1R, pick_vol
    decide_rbs2 = True
    find_po = pak_fam.index(str(name_capL1))  # find the position of the selected capsid in the name list
    trip = find_po  # equivalent matching pos in other lists items
    if trip < 2:
        original_voly = vol_outerradR[trip]
        find_ogvol.append(original_voly)
        genometypeR = gen_vol
        gen_spR = genometypeR * (2 / 3)  # convert to the volume of a sphere
        gen_icosR = gen_spR * (0.983631649)  # convert to a icosahedron
        # convert this to a capsid with a new overall outer radius(this will be different for each capsid - scaling value)
        gen_changR = gen_icosR * scalR[
            trip]  # for the original conversion between inner and outer capsid volumes for each virus
        volume_desiredR = gen_changR  # the volume of our gene will determines the minimum required/target volume of our capsid
        vol_outerdes.append(volume_desiredR)
        pi = 3.14159265359
        # find the area of the triangular face
        rmR = 0.8090169944 * ((volume_desiredR / 2.181694991) ** (1. / 3.))  # find a radius for the target capsid
        side_lengthR = rmR / (0.80901699)
        total_areaR = (8.6602540378) * ((side_lengthR) ** 2)
        area_faceR = total_areaR / 20  # single face
        store_faceR.append(area_faceR)
        print("The area of a single triangular face from your required icosahedral capsid is: " + str(
            int(area_faceR)) + " A^2.")
        print(" Which had a volume of: " + str(gen_vol) + " A^3")
    elif trip >= 2:
        original_voly = vol_outerradR[trip]
        find_ogvol.append(original_voly)
        genometypeR1 = gen_volss
        gen_spR = genometypeR1 * (2 / 3)  # convert to the volume of a sphere
        gen_icosR = gen_spR * (0.983631649)  # convert to a icosahedron
        # convert this to a capsid with a new overall outer radius(this will be different for each capsid - scaling value)
        gen_changR = gen_icosR * scalR[
            trip]  # for the original conversion between inner and outer capsid volumes for each virus
        volume_desiredR = gen_changR  # the volume of our gene will determines the minimum required/target volume of our capsid
        vol_outerdes.append(volume_desiredR)
        pi = 3.14159265359
        # find the area of the triangular face
        rmR = 0.8090169944 * ((volume_desiredR / 2.181694991) ** (1. / 3.))  # find a radius for the target capsid
        side_lengthR = rmR / (0.80901699)
        total_areaR = (8.6602540378) * ((side_lengthR) ** 2)
        area_faceR = total_areaR / 20  # single face
        store_faceR.append(area_faceR)
        print("The area of a single triangular face from your required icosahedral capsid is: " + str(
            int(area_faceR)) + " A^2.")
        print(" Which had a volume of: " + str(gen_volss) + " A^3")

    # find the volume of the original capsid
    volume_ic1 = low_redVol  # volume of og capsid from the outer radius value - applicable to the change in tri facet size
    seeker = find_PC_again.index(volume_ic1)
    volume_ic1R = vol_outerradR[int(seeker)]
    rm_new1R = 0.8090169944 * ((volume_ic1R / 2.181694991) ** (1. / 3.))
    side_length1R = rm_new1R / (0.80901699)
    total_area1R = (8.6602540378) * ((side_length1R) ** 2)
    area_face1R = total_area1R / 20  # area of the original side face
    print("The area of a single triangular face of the " + str(name_capL1) + "viral vector is: "
          + str(int(area_face1R)) + " A^2." + " Which had a volume of: " + str(int(low_redVol)) + " A")
    # find the new T value
    store_newTR = []
    store_newprot = []
    store_enlg = []
    oldT = round(int(T_famR[trip]))
    oldprot = round((int(sub_ogR[trip])))
    for x in store_faceR:
        enlargeR = int(x) / area_face1R
        store_enlg.append(enlargeR)
        newTR = round(int(oldT) * enlargeR)
        store_newTR.append(newTR)  # hold all the new t values
        num_protein1R = 60 * (newTR)
        store_newprot.append(num_protein1R)
        # this will store the new number of protein required
    enlargRatio = store_enlg[0]
    new_amountofprotein = store_newprot[0]
    pick_vol = vol_outerdes[0]
    more_protein = int(new_amountofprotein) - int(oldprot)
    for ttr in store_newTR:
        if enlargRatio < 1:
            print("\nThe enlargement ratio between your required capsid and selected capsid is: " + str(enlargRatio) +
                  "\nThis means the new capsid has the same triangulation number of the original capsid which was " + str(oldT))
            print("\n*From this result, we can confirm that this viral capsid is most appropriate choice for the delivery of your gene and does not require further readjustments.")
            print("\n\n\t\t\t\t\t\t\t_Summary of Viral Capsid Optimization_\n")
            titles = ['Selected Capsid', 'Capsid Volume(A)', 'T Value(A)',
                      'No. Subunits']
            data = [name_capL1, volume_ic1R, oldT, oldprot]
            print("\t" + "|  " + "    |   ".join(titles) + " |")
            print("\t" + "-" * 87)
            print("\t   " + "\t\t\t\t\t".join(map(str, data)))
        elif enlargRatio >= 1:
            print("\nThe enlargement ratio between your required capsid and selected capsid is: " + str(
                enlargRatio) + "\nThis means the new capsid has a triangulation number of " + str(ttr))
            print("To deliver your gene, the original " + str(
                name_capL1) + "viral vector must be increased in size by a factor of " + str(enlargRatio))
            print("In order to acheive this, the number of proteins that make up this viral capsid will need to increase by " + str(round(more_protein)))
            print("\n\n\t\t\t\t\t\t\t_Summary of Viral Capsid Optimization_\n")
            titles = ['Selected Capsid', ' Original Capsid Volume(A)', 'Original T Value(A)',
                      'Original No. Subunits', 'New Packaging Volume(A)', 'New T Value(A)', 'New No. Subunits']
            data = [name_capL1, volume_ic1R, oldT, oldprot, pick_vol, ttr, new_amountofprotein]
            print("|" + "  | ".join(titles) + "|")
            print("-" * 180)
            print("\t\t\t\t\t ".join(map(str, data)))

def TVALreadjustmentsall():
    global vol_outerdesB, vol_outerradRB, store_newprotRB
    vol_outerdesB = []
    store_faceRB = []
    store_face2RB = []
    for tripB in range(0, 4):
        if tripB < 2:
            gen_spRB = gen_vol * (2 / 3)  # convert to the volume of a sphere
            gen_icosRB = gen_spRB * (0.983631649)  # convert to a icosahedron
            # convert this to a capsid with a new overall outer radius(this will be different for each capsid - scaling value)
            gen_changRB = gen_icosRB * scalRB[
                tripB]  # for the original conversion between inner and outer capsid volumes for each virus
            volume_desiredRB = gen_changRB  # the volume of our gene will determines the minimum required/target volume of our capsid
            vol_outerdesB.append(volume_desiredRB)
            pi = 3.14159265359
            # find the area of the triangular face
            rmRB = 0.8090169944 * ((volume_desiredRB / 2.181694991) ** (1. / 3.))  # find a radius for the target capsid
            side_lengthRB = rmRB / (0.80901699)
            total_areaRB = (8.6602540378) * ((side_lengthRB) ** 2)
            area_faceRB = total_areaRB / 20  # single face
            store_faceRB.append(area_faceRB)


        elif tripB >= 2:
            gen_spRB = gen_volss * (2 / 3)  # convert to the volume of a sphere
            gen_icosRB = gen_spRB * (0.983631649)  # convert to a icosahedron
            # convert this to a capsid with a new overall outer radius(this will be different for each capsid - scaling value)
            gen_changRB = gen_icosRB * scalRB[
                tripB]  # for the original conversion between inner and outer capsid volumes for each virus
            volume_desiredRB = gen_changRB  # the volume of our gene will determines the minimum required/target volume of our capsid
            vol_outerdesB.append(volume_desiredRB)
            pi = 3.14159265359
            # find the area of the triangular face
            rmRB = 0.8090169944 * ((volume_desiredRB / 2.181694991) ** (1. / 3.))  # find a radius for the target capsid
            side_lengthRB = rmRB / (0.80901699)
            total_areaRB = (8.6602540378) * ((side_lengthRB) ** 2)
            area_faceRB = total_areaRB / 20  # single face
            store_faceRB.append(area_faceRB)


        volume_ic1RB = vol_outerradRB[
            tripB]  # volume of og capsid from the outer radius value - applicable to the change in tri facet size
        rm_new1RB = 0.8090169944 * ((volume_ic1RB / 2.181694991) ** (1. / 3.))
        side_length1RB = rm_new1RB / (0.80901699)
        total_area1RB = (8.6602540378) * ((side_length1RB) ** 2)
        area_face1RB = total_area1RB / 20  # area of the original side face
        store_face2RB.append(area_face1RB)
    trackRB = -1
    store_newTRB = []
    store_newprotRB = []

    for x, y in zip(store_faceRB, store_face2RB):
        trackRB += 1
        enlargeRB = int(x) / int(y)
        newTRB = round(int(T_famRB[trackRB]) * enlargeRB)
        if int(newTRB) <= int(T_famRB[trackRB]):
            newTRB = int((T_famRB[trackRB]))
            store_newTRB.append(newTRB)  # hold all the new t values
            num_protein1RB = 60 * (newTRB)
            store_newprotRB.append(num_protein1RB)  # this will store the new number of protein required
        elif int(newTRB) > int(T_famRB[trackRB]):
            store_newTRB.append(newTRB)  # hold all the new t values
            num_protein1RB = 60 * (newTRB)
            store_newprotRB.append(num_protein1RB)  # this will store the new number of protein required

    print("\n" * 3 + "\t" * 15 + "_FINAL OVERALL SUMMARY TABLE OF VIRAL CAPSID OPTIMIZATION _\n")
    titles = ['Viral Capsid', 'Packaging Volume(A)', 'T Value(A)',
              'No. Subunits', 'New Volume(A)', 'New T Value(A)', 'New No. Subunits']
    data11 = [titles] + list(
        zip(pak_famB, vol_outerradRB, T_famRB, sub_ogRB, vol_outerdesB, store_newTRB, store_newprotRB))
    for i, d in enumerate(data11):
        line = ''.join(str(x).ljust(24) for x in d)
        print("\t" + line)
        if i == 0:
            print('\t' + '-' * len(line) + '\n')

def capsidstabwithAVAL():
    tirpos = pak_fam.index(str(name_cap))
    pos_stab = tirpos
    global capsid_stab_og
    if pos_stab == 1:
        capsid_stab_og = 0
        print("There is no data available for capsid stability of the Herpesviridae viral capsid.\n")
        # print capsid stability values for other viral vectors
        titles = ['Family', 'Outer Diameter', 'Sum Association Energy','Capsid Stability']
        data11 = [titles] + list(zip(pk_lis, outer_diam, energies_vir, stability_vir))
        for i, d in enumerate(data11):
            line = ''.join(str(x).ljust(24) for x in d)
            print("\t" + line)
            if i == 0:
                print('\t' + '-' * len(line) + '\n')
    elif pos_stab == 0:
        capsid_stab_og = stability_vir[0]
        print("The " + str(name_cap) + " viral capsid has an overall volume of " + str(original_voly) + " A.")
        print("The overall capsid stability was calcualted from the summation of association energies at each unique interface of protein subunits, and noramalied by its respective diameter")
        print("THe capsid stability of the " + str(name_cap) + "viral vector is equal to:" + str(capsid_stab_og) + "kJ/A")
    elif pos_stab == 2:
        capsid_stab_og = stability_vir[1]
        print("The " + str(name_cap) + " viral capsid has an overall volume of " + str(original_voly) + " A.")
        print("The overall capsid stability was calcualted from the summation of association energies at each unique interface of protein subunits, and noramalied by its respective diameter")
        print("THe capsid stability of the " + str(name_cap) + "viral vector is equal to:" + str(capsid_stab_og) + "kJ/A")
    elif pos_stab == 3:
        capsid_stab_og = stability_vir[2]
        print("The " + str(name_cap) + " viral capsid has an overall volume of " + str(original_voly) + " A.")
        print("The overall capsid stability was calcualted from the summation of association energies at each unique interface of protein subunits, and noramalied by its respective diameter")
        print("THe capsid stability of the " + str(name_cap) + "viral vector is equal to:" + str(capsid_stab_og) + "kJ/A")

def capsidstabwithNOOAVAL():
    tirpos2 = pak_fam.index(str(name_capL1))
    pos_stab = tirpos2
    global capsid_stab_og, new_capstability
    if pos_stab == 1:
        capsid_stab_og = 0
        print("There is no data available for capsid stability of the Herpesviridae viral capsid. ")
        titles = ['Family', 'Outer Diameter', 'Sum Association Energy',
                  'Capsid Stability']
        data11 = [titles] + list(zip(pk_lis, outer_diam, energies_vir, stability_vir))
        for i, d in enumerate(data11):
            line = ''.join(str(x).ljust(24) for x in d)
            print("\t" + line)
            if i == 0:
                print('\t' + '-' * len(line) + '\n')
    elif pos_stab == 0:
        capsid_stab_og = stability_vir[0]
        print("The " + str(name_capL1) + " viral capsid has an overall volume of " + str(original_voly) + " A.")
        print(
            "Capsid stability was calculated from the summation of association energies at each unique interface of protein subunits, and normalized by its respective diameter")
        print("\nThe capsid stability of the " + str(name_capL1) + "viral vector is equal to: " + str(
            capsid_stab_og) + " kJ/A")
        # find the sum of the association energies - selct capsid
        pos_eng = stability_vir.index(capsid_stab_og)
        cap_engdec = energies_vir[pos_eng]
        # find hte new diameter and calculate the new capsid stability
        cap_rad_theo = rmR
        new_diam = int(cap_rad_theo) * 2
        new_capstability = int(cap_engdec) / new_diam
        print("The volume of the capsid required for your gene is: " + str(
            volume_desiredR) + " A." + " This has a capsid stability of: " + str(new_capstability) + " kJ/A")
        # compare with the old stability value
        if new_capstability > capsid_stab_og:
            print("Original Capsid Stability (kJ/A) : " + str(capsid_stab_og))
            print("New Capsid Stability (kJ/A) : " + str(new_capstability))
            print("The original capsid for the " + str(
                name_capL1) + "viral vector is more stable than the new capsid constructed for your selected gene\n\n")
        elif new_capstability < capsid_stab_og:
            print("Original Capsid Stability: " + str(capsid_stab_og))
            print("New Capsid Stability: " + str(new_capstability))
            print("The new capsid constructed specifically for your gene from the " + str(
                name_capL1) + "viral vector, is more stable than the original capsid\n\n")
        titles = ['Family', 'Outer Diameter', 'Sum Association Energy',
                  'Capsid Stability']
        data11 = [titles] + list(zip(pk_lis, outer_diam, energies_vir, stability_vir))
        for i, d in enumerate(data11):
            line = ''.join(str(x).ljust(24) for x in d)
            print("\t" + line)
            if i == 0:
                print('\t' + '-' * len(line) + '\n')

    elif pos_stab == 2:
        capsid_stab_og = stability_vir[1]
        print("The " + str(name_capL1) + " viral capsid has an overall volume of " + str(original_voly) + " A.")
        print(
            "Capsid stability was calculated from the summation of association energies at each unique interface of protein subunits, and normalized by its respective diameter")
        print("\nThe capsid stability of the " + str(name_capL1) + "viral vector is equal to: " + str(
            capsid_stab_og) + " kJ/A")
        # find the sum of the association energies - selct capsid
        pos_eng = stability_vir.index(capsid_stab_og)
        cap_engdec = energies_vir[pos_eng]
        # find hte new diameter and calculate the new capsid stability
        cap_rad_theo = rmR
        new_diam = int(cap_rad_theo) * 2
        new_capstability = int(cap_engdec) / new_diam
        print("The volume of the capsid required for your gene is: " + str(
            volume_desiredR) + " A." + " This has a capsid stability of: " + str(new_capstability) + " kJ/A")
        # compare with the old stability value
        if new_capstability > capsid_stab_og:
            print("Original Capsid Stability (kJ/A) : " + str(capsid_stab_og))
            print("New Capsid Stability (kJ/A) : " + str(new_capstability))
            print("The original capsid for the " + str(
                name_capL1) + "viral vector is more stable than the new capsid constructed for your selected gene\n\n")
        elif new_capstability < capsid_stab_og:
            print("Original Capsid Stability: " + str(capsid_stab_og))
            print("New Capsid Stability: " + str(new_capstability))
            print("The new capsid constructed specifically for your gene from the " + str(
                name_capL1) + "viral vector, is more stable than the original capsid\n\n")
        titles = ['Family', 'Outer Diameter', 'Sum Association Energy',
                  'Capsid Stability']
        data11 = [titles] + list(zip(pk_lis, outer_diam, energies_vir, stability_vir))
        for i, d in enumerate(data11):
            line = ''.join(str(x).ljust(24) for x in d)
            print("\t" + line)
            if i == 0:
                print('\t' + '-' * len(line) + '\n')
    elif pos_stab == 3:
        capsid_stab_og = stability_vir[2]
        print("The " + str(name_capL1) + " viral capsid has an overall volume of " + str(original_voly) + " A.")
        print(
            "Capsid stability was calculated from the summation of association energies at each unique interface of protein subunits, and normalized by its respective diameter")
        print("\nThe capsid stability of the " + str(name_capL1) + "viral vector is equal to: " + str(
            capsid_stab_og) + " kJ/A")
        # find the sum of the association energies - selct capsid
        pos_eng = stability_vir.index(capsid_stab_og)
        cap_engdec = energies_vir[pos_eng]
        # find hte new diameter and calculate the new capsid stability
        cap_rad_theo = rmR
        new_diam = int(cap_rad_theo) * 2
        new_capstability = int(cap_engdec) / new_diam
        print("The volume of the capsid required for your gene is: " + str(
            volume_desiredR) + " A." + " This has a capsid stability of: " + str(new_capstability) + " kJ/A")
        # compare with the old stability value
        if new_capstability > capsid_stab_og:
            print("Original Capsid Stability (kJ/A) : " + str(capsid_stab_og))
            print("New Capsid Stability (kJ/A) : " + str(new_capstability))
            print("The original capsid for the " + str(
                name_capL1) + "viral vector is more stable than the new capsid constructed for your selected gene\n\n")
        elif new_capstability < capsid_stab_og:
            print("Original Capsid Stability: " + str(capsid_stab_og))
            print("New Capsid Stability: " + str(new_capstability))
            print("The new capsid constructed specifically for your gene from the " + str(
                name_capL1) + "viral vector, is more stable than the original capsid\n\n")
        titles = ['Family', 'Outer Diameter', 'Sum Association Energy',
                  'Capsid Stability']
        data11 = [titles] + list(zip(pk_lis, outer_diam, energies_vir, stability_vir))
        for i, d in enumerate(data11):
            line = ''.join(str(x).ljust(24) for x in d)
            print("\t" + line)
            if i == 0:
                print('\t' + '-' * len(line) + '\n')
def FINALTABLE1():
    print("\n" * 3 + "\t" * 18 + "_SUMMARY OF VIRAL CAPSID OPTIMIZATION_\n")
    titles = ['Selected Capsid','Packing Volume', 'T value', 'Subunits Proteins', 'Capsid Stability']
    data = [name_cap, round(volume_ic1R), oldT, oldprot, capsid_stab_og]
    print("\t" * 12 + "|" + "  | ".join(titles) + "|")
    print("\t" * 12 +"-" * 90)
    print('\t '*12 + "\t\t\t".join(map(str, data)))
def FINALTABLE2():
    print("\n" * 3 + "\t" * 18 + "_SUMMARY OF VIRAL CAPSID OPTIMIZATION_\n")
    titles = ['Selected Capsid','Packing Volume', 'T value', 'Subunits Proteins', 'Capsid Stability']
    data = [name_capL1, round(pick_vol), ttr, new_amountofprotein, capsid_stab_og]
    print("\t" * 12 + "|" + "  | ".join(titles) + "|")
    print("\t" * 12 + "-" * 90)
    print('\t ' * 12 + "\t\t\t".join(map(str, data)))

def therapeuticsTABLE():
    # SUMMARY TABLE BOTTOM
    pack_cap = [7.5, 7.5, ">30", ">30", 4.5, 4.5, 4.5, 8, 8, 7.5]
    drugs = [["Gendicine", "Oncorine"], ["Imlygic", "Zalmoxis", "Neovaculgen"], ["Luxturna", "Zolgensma", "Spinraza"],
             ["Kymriah", "Yescarta", "Zalmoxis"], ["None"]]

    # release drugs into one list
    drug_lis = []
    for row in drugs:
        for elem in row:
            drug_lis.append(elem)

    # print only average capsid volumes for parvoviridae
    avg_parvov = []
    for parvolE in parvov_pos:
        parv_vol = capsid_final[parvolE]
        avg_parvov.append(parv_vol)
        avg_P = sum(avg_parvov) / len(avg_parvov)  # this is the new value for the parvov virus capsid volume

    # create new list for the table, that will only include the one avg value for the
    capvol_RT = [145124724.904, 145124724.904, 384253856.516159, 384253856.516159, avg_P, avg_P, avg_P, 817283.234436,
                 817283.234436, 77072667.4378]
    fam_listRT = ['Adenoviridae', 'Adenoviridae', 'Herpesviridae', 'Herpesviridae', 'Parvoviridae', 'Parvoviridae',
                  'Parvoviridae', 'Retroviridae', 'Retroviridae', 'Togaviridae']
    vec = ["rAD-p53", "H103", "pCMV-vegf165", "HIV-1", "AAV2", "AAV9", "AAV9", "y_RV", "MLV(?)"]
    pdb_num = ['6CGU', '6CGU', '1CMV', 'IDLO', '1LP3', '3UX1', '3UX1', '?', '6HWW']
    url_L = ['http://www.rcsb.org/structure/6CGV', 'http://www.rcsb.org/structure/6CGV',
             'http://www.rcsb.org/structure/1CMV',
             'http://www.rcsb.org/structure/1DLO', 'https://www.rcsb.org/structure/1lp3',
             'https://www.rcsb.org/structure/3UX1',
             'https://www.rcsb.org/structure/3UX1', '?', 'http://www.rcsb.org/structure/6HWW']

    # TABLE bottom
    print("\n" * 3 + "\t" * 15 + "_COMPLETE TABLE OF VIRAL VECTORS FOR GENE THERAPY_\n")
    titles = ['Family', 'Capsid Volume', 'Packing Capacity(kb)', 'Approved Therapy', 'Viral Vector', 'PDB Number',
              'Link']
    data = [titles] + list(zip(fam_listRT, capvol_RT, pack_cap, drug_lis, vec, pdb_num, url_L))

    for i, d in enumerate(data):
        line = ''.join(str(x).ljust(22) for x in d)
        print("\t" + line)
        if i == 0:
            print('\t' + '-' * len(line) + '\n')
    print("\n**These include viruses used in approved therapies and ones that are being researched in pre-clinical and clinical trials**\n\n\n")

def sortandfindPCvals2():
    global pp1, pp2, pp1_V, pp2_V
    count_pack = -1  # this will record the names of positions in both conditions
    pp1 = []  # will hold family names of the available capsids
    pp2 = []  # will hold the family names of the unavailable capsids ---- both are based on packing volumes(cylindrical)
    pp1_V = []
    pp2_V = []  # this will hold the volumes of the unavailable capsids
    pak_fam = ['Adenoviridae', 'Herpesviridae', 'Parvoviridae', 'Retroviridae']
    T_fam = [25, 16, 1, 1]  # holds the T values of the approved capsids
    counter_packingV = 0  # this used for output statements later
    case1 = False
    case2 = False
    re1 = 0  # this will record the positions of the ones that are approved
    re1store = []  # this will store all the available positions
    # Unpack volumes of the calculated volumes within a loop
    for find1 in capsVlis[0:2]:  # between only the ds capsid volumes
        if int(gen_vol) < int(
                find1):  # THIS CAN BE USED FOR THE GENE - APPROVED PACKING CAPACITY - this compares a ds gene to capsV_lis - should change according to same genome type
            counter_packingV += 1  # this is used later for the output list
            count_pack += 1  # record the position inside the list ---use this to find the family namesw
            new_pack2 = pak_fam[count_pack]  # - use the position to find the name of the capsid
            pp1.append(new_pack2)  # collect all the names of the capsids larger than cargo ---can be used
            pp1_V.append(find1)  # this will collect only available volumes (calc of pack)
            case1 = True
            re1 += 1
            re1store.append(
                re1)  # hold onto all the positions that are added -- this order will be printed in the list below
        if int(gen_vol) > int(find1):  # THIS CAN NOT BE USED FOR THE GENE ---APPROVED PACKING CAPACITY
            count_pack += 1
            new_pack = pak_fam[count_pack]
            pp2.append(new_pack)  # collect all the capsids smaller than cargo -- can not be used
            pp2_V.append(find1)
            case2 = True
        if count_pack == 1:  # when the loop is finished evaluating all the ds genome type then go to the next two ss genome types
            for find2 in capsVlis[2:5]:
                if int(gen_volss) < find2:
                    counter_packingV += 1
                    count_pack += 1
                    new_pack3 = pak_fam[count_pack]
                    pp1.append(new_pack3)  # collect the names of the ss capsids
                    pp1_V.append(find2)  # collect the volumes of ss capsid volumes
                    case1 = True
                    re1 += 1
                    re1store.append(re1)
                if int(gen_volss) > find2:
                    count_pack += 1
                    new_pack4 = pak_fam[count_pack]
                    pp2.append(new_pack4)  # add the capsid name to the pp2 list
                    pp2_V.append(find2)  # this will collect the capsids volumes that are unavailable
                    case2 = True
def havePCatleastone2():
    global name_cap, smallest_volPC
    # within the pack table, find the smallest capsid ---- there must be at least one entry inside this table - this will be the selected capsid
    smallest_volPC = min(pp1_V)  # this will hold the smallest volume
    name_cap_pos = pp1_V.index(smallest_volPC)
    name_cap = pp1[name_cap_pos]

def noPCone2():
    global name_capL1, low_redVol
    adj_pcVf = []
    adj_pcVi = []
    for vee in capsVlis:
        if int(vee) == 8010000:
            diff_pcVee6 = gen_vol - int(vee)  # only ds nucleic acid
            add_pcVee6 = int(vee) + diff_pcVee6
            adj_pcVi.append(diff_pcVee6)
            adj_pcVf.append(add_pcVee6)
        if int(vee) == 32040000:
            diff_pcVee7 = gen_vol - int(vee)  # only ds nucleic acid
            add_pcVee7 = int(vee) + diff_pcVee7
            adj_pcVi.append(diff_pcVee7)
            adj_pcVf.append(add_pcVee7)
        if int(vee) == 2403000:
            diff_pcVee8 = gen_volss - int(vee)
            add_pcVee8 = int(vee) + diff_pcVee8
            adj_pcVi.append(diff_pcVee8)
            adj_pcVf.append(add_pcVee8)
        if int(vee) == 4272000:
            diff_pcVee9 = gen_volss - int(vee)
            add_pcVee9 = int(vee) + diff_pcVee9
            adj_pcVi.append(diff_pcVee9)
            adj_pcVf.append(add_pcVee9)
    # find the one requires the smallest readjustment
    low_red = min(adj_pcVi)
    name_capL2 = adj_pcVi.index(low_red)  # holds position in a list of 4
    name_capL1 = pp2[name_capL2]
    low_redVol = capsVlis[name_capL2]


def availablePTCVAL2():
    vol_outerdes = []
    find_ogvol = []
    store_faceR = []
    gen_typer = []
    global decide_rbs1, trip, original_voly, rmR, volume_desiredR
    decide_rbs1 = True
    find_po = pak_fam.index(str(name_cap))  # find the position of the selected capsid in the name list
    trip = find_po  # equivalent matching pos in other lists items
    if trip < 2:
        original_voly = vol_outerradR[trip]
        find_ogvol.append(original_voly)
        genometypeR = gen_vol
        gen_typer.append(genometypeR)
        gen_spR = genometypeR * (2 / 3)  # convert to the volume of a sphere
        gen_icosR = gen_spR * (0.983631649)  # convert to a icosahedron
        # convert this to a capsid with a new overall outer radius(this will be different for each capsid - scaling value)
        gen_changR = gen_icosR * scalR[
            trip]  # for the original conversion between inner and outer capsid volumes for each virus
        volume_desiredR = gen_changR  # the volume of our gene will determines the minimum required/target volume of our capsid
        vol_outerdes.append(volume_desiredR)
        pi = 3.14159265359
        # find the area of the triangular face
        rmR = 0.8090169944 * ((volume_desiredR / 2.181694991) ** (1. / 3.))  # find a radius for the target capsid
        side_lengthR = rmR / (0.80901699)
        total_areaR = (8.6602540378) * ((side_lengthR) ** 2)
        area_faceR = total_areaR / 20  # single face
        store_faceR.append(area_faceR)

    elif trip >= 2:
        original_voly = vol_outerradR[trip]
        find_ogvol.append(original_voly)
        genometypeR1 = gen_volss
        gen_typer.append(genometypeR1)
        gen_spR = genometypeR1 * (2 / 3)  # convert to the volume of a sphere
        gen_icosR = gen_spR * (0.983631649)  # convert to a icosahedron
        # convert this to a capsid with a new overall outer radius(this will be different for each capsid - scaling value)
        gen_changR = gen_icosR * scalR[
            trip]  # for the original conversion between inner and outer capsid volumes for each virus
        volume_desiredR = gen_changR  # the volume of our gene will determines the minimum required/target volume of our capsid
        vol_outerdes.append(volume_desiredR)
        pi = 3.14159265359
        # find the area of the triangular face
        rmR = 0.8090169944 * ((volume_desiredR / 2.181694991) ** (1. / 3.))  # find a radius for the target capsid
        side_lengthR = rmR / (0.80901699)
        total_areaR = (8.6602540378) * ((side_lengthR) ** 2)
        area_faceR = total_areaR / 20  # single face
        store_faceR.append(area_faceR)


    # find the volume of the original capsid
    volume_ic1 = smallest_volPC
    seeker = find_PC_again.index(volume_ic1)
    volume_ic1R = vol_outerradR[int(seeker)]
    rm_new1R = 0.8090169944 * ((volume_ic1R / 2.181694991) ** (1. / 3.))
    side_length1R = rm_new1R / (0.80901699)
    total_area1R = (8.6602540378) * ((side_length1R) ** 2)
    area_face1R = total_area1R / 20  # area of the original side face
    # find the new T value
    store_newTR = []
    store_newprot = []
    store_enlg = []
    oldT = round(int(T_famR[trip]))
    oldprot = round((int(sub_ogR[trip])))
    for x in store_faceR:
        enlargeR = int(x) / area_face1R
        store_enlg.append(enlargeR)
        newTR = round(int(oldT) * enlargeR)
        store_newTR.append(newTR)  # hold all the new t values
        num_protein1R = 60 * (newTR)
        store_newprot.append(num_protein1R)
        # this will store the new number of protein required
    enlargRatio = store_enlg[0]
    new_amountofprotein = store_newprot[0]
    more_protein = int(new_amountofprotein) - int(oldprot)

def NOavailablePCTVAL2():
    vol_outerdes = []
    find_ogvol = []
    store_faceR = []
    gen_typer = []
    global decide_rbs2, trip, original_voly, rmR, volume_desiredR
    decide_rbs2 = True
    find_po = pak_fam.index(str(name_capL1))  # find the position of the selected capsid in the name list
    trip = find_po  # equivalent matching pos in other lists items
    if trip < 2:
        original_voly = vol_outerradR[trip]
        find_ogvol.append(original_voly)
        genometypeR = gen_vol
        gen_spR = genometypeR * (2 / 3)  # convert to the volume of a sphere
        gen_icosR = gen_spR * (0.983631649)  # convert to a icosahedron
        # convert this to a capsid with a new overall outer radius(this will be different for each capsid - scaling value)
        gen_changR = gen_icosR * scalR[
            trip]  # for the original conversion between inner and outer capsid volumes for each virus
        volume_desiredR = gen_changR  # the volume of our gene will determines the minimum required/target volume of our capsid
        vol_outerdes.append(volume_desiredR)
        pi = 3.14159265359
        # find the area of the triangular face
        rmR = 0.8090169944 * ((volume_desiredR / 2.181694991) ** (1. / 3.))  # find a radius for the target capsid
        side_lengthR = rmR / (0.80901699)
        total_areaR = (8.6602540378) * ((side_lengthR) ** 2)
        area_faceR = total_areaR / 20  # single face
        store_faceR.append(area_faceR)

    elif trip >= 2:
        original_voly = vol_outerradR[trip]
        find_ogvol.append(original_voly)
        genometypeR1 = gen_volss
        gen_spR = genometypeR1 * (2 / 3)  # convert to the volume of a sphere
        gen_icosR = gen_spR * (0.983631649)  # convert to a icosahedron
        # convert this to a capsid with a new overall outer radius(this will be different for each capsid - scaling value)
        gen_changR = gen_icosR * scalR[
            trip]  # for the original conversion between inner and outer capsid volumes for each virus
        volume_desiredR = gen_changR  # the volume of our gene will determines the minimum required/target volume of our capsid
        vol_outerdes.append(volume_desiredR)
        pi = 3.14159265359
        # find the area of the triangular face
        rmR = 0.8090169944 * ((volume_desiredR / 2.181694991) ** (1. / 3.))  # find a radius for the target capsid
        side_lengthR = rmR / (0.80901699)
        total_areaR = (8.6602540378) * ((side_lengthR) ** 2)
        area_faceR = total_areaR / 20  # single face
        store_faceR.append(area_faceR)

    # find the volume of the original capsid
    volume_ic1 = low_redVol  # volume of og capsid from the outer radius value - applicable to the change in tri facet size
    seeker = find_PC_again.index(volume_ic1)
    volume_ic1R = vol_outerradR[int(seeker)]
    rm_new1R = 0.8090169944 * ((volume_ic1R / 2.181694991) ** (1. / 3.))
    side_length1R = rm_new1R / (0.80901699)
    total_area1R = (8.6602540378) * ((side_length1R) ** 2)
    area_face1R = total_area1R / 20  # area of the original side face

    # find the new T value
    store_newTR = []
    store_newprot = []
    store_enlg = []
    oldT = round(int(T_famR[trip]))
    oldprot = round((int(sub_ogR[trip])))
    for x in store_faceR:
        enlargeR = int(x) / area_face1R
        store_enlg.append(enlargeR)
        newTR = round(int(oldT) * enlargeR)
        store_newTR.append(newTR)  # hold all the new t values
        num_protein1R = 60 * (newTR)
        store_newprot.append(num_protein1R)
        # this will store the new number of protein required
    enlargRatio = store_enlg[0]
    new_amountofprotein = store_newprot[0]
    more_protein = int(new_amountofprotein) - int(oldprot)

def TVALreadjustmentsall2():
    global store_newprotRB, vol_outerdesB, store_newTRB
    vol_outerdesB = []
    store_faceRB = []
    store_face2RB = []
    for tripB in range(0, 4):
        if tripB < 2:
            gen_spRB = gen_vol * (2 / 3)  # convert to the volume of a sphere
            gen_icosRB = gen_spRB * (0.983631649)  # convert to a icosahedron
            # convert this to a capsid with a new overall outer radius(this will be different for each capsid - scaling value)
            gen_changRB = gen_icosRB * scalRB[
                tripB]  # for the original conversion between inner and outer capsid volumes for each virus
            volume_desiredRB = gen_changRB  # the volume of our gene will determines the minimum required/target volume of our capsid
            vol_outerdesB.append(volume_desiredRB)
            pi = 3.14159265359
            # find the area of the triangular face
            rmRB = 0.8090169944 * ((volume_desiredRB / 2.181694991) ** (1. / 3.))  # find a radius for the target capsid
            side_lengthRB = rmRB / (0.80901699)
            total_areaRB = (8.6602540378) * ((side_lengthRB) ** 2)
            area_faceRB = total_areaRB / 20  # single face
            store_faceRB.append(area_faceRB)

        elif tripB >= 2:
            gen_spRB = gen_volss * (2 / 3)  # convert to the volume of a sphere
            gen_icosRB = gen_spRB * (0.983631649)  # convert to a icosahedron
            # convert this to a capsid with a new overall outer radius(this will be different for each capsid - scaling value)
            gen_changRB = gen_icosRB * scalRB[
                tripB]  # for the original conversion between inner and outer capsid volumes for each virus
            volume_desiredRB = gen_changRB  # the volume of our gene will determines the minimum required/target volume of our capsid
            vol_outerdesB.append(volume_desiredRB)
            pi = 3.14159265359
            # find the area of the triangular face
            rmRB = 0.8090169944 * ((volume_desiredRB / 2.181694991) ** (1. / 3.))  # find a radius for the target capsid
            side_lengthRB = rmRB / (0.80901699)
            total_areaRB = (8.6602540378) * ((side_lengthRB) ** 2)
            area_faceRB = total_areaRB / 20  # single face
            store_faceRB.append(area_faceRB)

        volume_ic1RB = vol_outerradRB[
            tripB]  # volume of og capsid from the outer radius value - applicable to the change in tri facet size
        rm_new1RB = 0.8090169944 * ((volume_ic1RB / 2.181694991) ** (1. / 3.))
        side_length1RB = rm_new1RB / (0.80901699)
        total_area1RB = (8.6602540378) * ((side_length1RB) ** 2)
        area_face1RB = total_area1RB / 20  # area of the original side face
        store_face2RB.append(area_face1RB)
    trackRB = -1
    store_newTRB = []
    store_newprotRB = []

    for x, y in zip(store_faceRB, store_face2RB):
        trackRB += 1
        enlargeRB = int(x) / int(y)
        newTRB = round(int(T_famRB[trackRB]) * enlargeRB)
        if int(newTRB) <= int(T_famRB[trackRB]):
            newTRB = int((T_famRB[trackRB]))
            store_newTRB.append(newTRB)  # hold all the new t values
            num_protein1RB = 60 * (newTRB)
            store_newprotRB.append(num_protein1RB)  # this will store the new number of protein required
        elif int(newTRB) > int(T_famRB[trackRB]):
            store_newTRB.append(newTRB)  # hold all the new t values
            num_protein1RB = 60 * (newTRB)
            store_newprotRB.append(num_protein1RB)  # this will store the new number of protein required

def calcRATIO():
    global all_liss
    import math
    import random
    li = []
    for _ in range(16):
        output = []
        n = int(store_newprotRB[2] / 20)
        cn = int(n) if ((int(n) % 2) == 0) else int(n) + 1
        newratA = random.randrange(0, cn + 1, 2)
        newratB = random.randrange(newratA, cn + 1, 2) - newratA
        newratC = cn - newratA - newratB
        output.extend(value for name, value in locals().items() if name.startswith('newrat'))
        li.append(output)
    fail = []
    passer = []
    counter = 0
    counter2 = 0
    for x in li:
        if 0 in x:
            fail.append(x)
        else:
            passer.append(x)
    all_liss = []
    for p in passer:
        liss = []
        res = math.gcd(*p[:2])  # get the gcd of first two numbers
        for x in p[2:]:  # now iterate over the list starting from the 3rd element
            res1 = math.gcd(res, x)
        for y in p:
            new_values = int(y / res1)
            liss.append(new_values)
        all_liss.append(liss)
    new_all_liss = []  # remove the repeated ratio lists
    for elem in all_liss:
        if elem not in new_all_liss:
            new_all_liss.append(elem)
    all_liss = new_all_liss
    print("\nThe original ratio of protein in the capsid: 1 : 1 : 10 ")
    print("\nNew potential ratio(s) of subunits proteins in the new capsid: ")
    for j in all_liss:
        print(' : '.join(map(repr, j)))
def ratioTABLE():
    print("\n" * 3 + "\t" * 18 + "_AAV VIRAL CAPSID OPTIMIZATION _\n")

    titles = ['Viral Capsid', 'Packaging Volume(A)', 'T Value(A)',
              'No. Subunits', 'New Volume(A)', 'New T Value(A)', 'New No. Subunits']
    data = ([pak_famB[2], vol_outerradRB[2], T_famRB[2], sub_ogRB[2], int(vol_outerdesB[2]), store_newTRB[2],
             store_newprotRB[2]])
    print("\t\t\t|\t " + "\t | \t ".join(titles) + "\t|")
    print(" " * 10 + "-" * 160)
    print("\t" * 4 + "\t\t\t\t".join(map(str, data)))

'''
plot graphs 
'''
def plotavailpie():
    import plotly.graph_objects as go

    labels = ['Available Capsids','Unavailable Capsids']
    values = [counter, 68-int(counter)]

    # Use `hole` to create a donut-like pie chart
    fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.3)])
    fig.update_layout(
        title_text="Proportion of Viral Capsid Availability",
       )
    fig.show()

def plotoptizgraph():
    import plotly.graph_objects as go
    li = vol_outerradRB
    y = []
    y.extend(li)
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=["Adenoviridae", 'Herpesviridae', 'Parvoviridae', 'Retroviridae'],
        y=y,
        name='Original Capsid Volume',
        marker_color='rgb(46, 48, 184)'
    ))
    ti = vol_outerdesB
    y = []
    y.extend(ti)
    fig.add_trace(go.Bar(
        x=["Adenoviridae", 'Herpesviridae', 'Parvoviridae', 'Retroviridae'],
        y=y,
        name='Optimized Capsid Volume',
        marker_color='rgb(105, 176, 219)'
    ))
    fig.update_layout(
        title='Optimized Packaging Volumes of Viral Capsids',
        xaxis_tickfont_size=14,
        yaxis=dict(
            title='Volume (A^3)',
            titlefont_size=16,
            tickfont_size=14,
        ))
    fig.show()
def conv_vol_diamNOCAPS():
    global diam_oldcap, diam_newcap
    pi = 3.14159265359
    oldcapvol = volume_ic1R
    diam_oldcap = round(((oldcapvol*6)/pi)**(1./3.))
    newcapvol = pick_vol
    diam_newcap = round(((newcapvol*6)/pi)**(1/3))
    return diam_oldcap, diam_newcap


def plot3dgraphNOCAPS():
    import plotly.graph_objects as go

    capsids = ['Optimised Capsid', 'Old Capsid']
    capsid_colors = ['rgb(135, 135, 125)', 'rgb(210, 50, 0)']

    capsid_diameter = [diam_newcap, diam_oldcap]


    fig = go.Figure(data=go.Scatter3d(
        x=[100, 100],
        y=[100, 100],
        z=[230.1, 10.0],
        text=capsids,
        mode='markers',
        marker=dict(
            sizemode='diameter',
            sizeref=10,  # info on sizeref: https://plot.ly/python/reference/#scatter-marker-sizeref
            size=capsid_diameter,
            color=capsid_colors,
        )
    ))

    fig.update_layout(width=600, height=600, title='Optimization or Viral Capsids')

    fig.show()

def conv_vol_diamAVCAPS(): #the capsids for this condition are the same, sine PK availability condtion is True
    global diam_oldcap, diam_newcap
    pi = 3.14159265359
    oldcapvol = volume_ic1R
    diam_oldcap = round(((oldcapvol*6)/pi)**(1./3.))
    newcapvol = volume_ic1R
    diam_newcap = round(((newcapvol*6)/pi)**(1/3))
    return diam_oldcap, diam_newcap

def plot3dgraphAVCAPS():
    import plotly.graph_objects as go
    capsids = ['Optimised Capsid', 'Old Capsid']
    capsid_colors = ['rgb(135, 135, 125)', 'rgb(210, 50, 0)']

    capsid_diameter = [diam_newcap, diam_oldcap]

    fig = go.Figure(data=go.Scatter3d(
        x=[100, 100],
        y=[100, 100],
        z=[230.1, 10.0],
        text=capsids,
        mode='markers',
        marker=dict(
            sizemode='diameter',
            sizeref=10,  # info on sizeref: https://plot.ly/python/reference/#scatter-marker-sizeref
            size=capsid_diameter,
            color=capsid_colors,
        )
    ))

    fig.update_layout(width=600, height=600, title='Optimization or Viral Capsids')

    fig.show()

def plotfinaltable1():
    import plotly.graph_objects as go
    fig = go.Figure(data=[go.Table(
        header=dict(values=['Selected Capsid', 'Packing Volume', ' T value', 'Subunit Proteins', 'Capsid Stability']),

        cells=dict(values=[[name_cap], [volume_ic1R], [oldT], [oldprot], [capsid_stab_og]]))
                          ])
    fig.show()

def plotfinaltable2():
    import plotly.graph_objects as go
    fig = go.Figure(data=[go.Table(
        header=dict(values=['Selected Capsid', 'Packing Volume', ' T value', 'Subunit Proteins', 'Capsid Stability']),

        cells=dict(values=[[name_capL1], [gen_vol], [ttr], [new_amountofprotein], [capsid_stab_og]]))
                          ])
    fig.show()

def summaryplotall():
    print(counter, diam_newcap, diam_oldcap)
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    capsids = ['Optimised Capsid', 'Original Capsid']
    capsid_colors = ['rgb(135, 135, 125)', 'rgb(210, 50, 0)']

    bar_y1 = store_newprotRB
    bar_y = []
    bar_y.extend(bar_y1)
    print(bar_y)


    labelsP = ['Available Capsids', 'Unavailable Capsids']
    valuesP = [counter, 68 - int(counter)]

    capsid_diameter = [diam_newcap, diam_oldcap]

    # Initialize figure with subplots
    fig = make_subplots(
        rows=3, cols=2,
        column_widths=[0.8, 0.6],
        row_heights=[0.01, 1.7, 1.0],
        specs=[[None, None], [{"type": "scatter3d"}, {"type": "bar"}],
               [None, {"type": "pie"}]])

    # Add locations bar chart
    fig.update_layout(
        title_text=" Viral Capsid Optimization", showlegend=False
    )
    fig.add_trace(go.Pie(labels=labelsP, values=valuesP, hole=.3),
                  row=3, col=2)

    fig.add_trace(go.Bar(x=['Adenoviridae', 'Herpesviridae', 'Parvoviridae', 'Retroviridae'], y=bar_y),
                  row=2, col=2)

    fig.add_trace(go.Scatter3d(x=[100, 100], y=[100, 100], z=[230.1, 10.0], text=capsids, mode='markers',
                               marker=dict(
                                   sizemode='diameter',
                                   sizeref=10,
                                   size=capsid_diameter,
                                   color=capsid_colors)),
                  row=2, col=1
                  )

    # Set theme, margin, and annotation in layout
    fig.update_layout(
        template="plotly_white",
        margin=dict(r=15, t=25, b=40, l=50),
        annotations=[
            go.layout.Annotation(
                text="Summary of Optimization Data: (1)Capsid Size Comparison Models (2) New Total Subunit Protein Levels (3) Total Capsid Availability (/68)",
                showarrow=False,
                xref="paper",
                yref="paper",
                x=0,
                y=0)
        ]
    )

    fig.show()



