'''
CapsidBuilder

Investigate the feasibility of novel viral capsid construction by evaluating the protein expression
levels from wet lab constructs using translation efficiency measurements. These measurements will be
used to select one of our constructs that may be to build novel theoretical viral capsids.

'''
from P1 import*
from P2 import*
import math
import plotly.graph_objects as go
from plotly.subplots import make_subplots
def script():
    table_presets()
    gene = ask_user()
    if gene not in gen_sizeR:
        print("There is no gene length that matches this in the table, please try again.")
        script()
    if gene == 0:  # if the user inputs 0, should restart script
        print("Your gene length must be greater than 0. Please try again.\n")
        script()

    # ds and ss genome volumes - cylindrical volumes
    gene_volumes = gene_vol()
    gen_vol = gene_volumes[0]
    gen_volss = gene_volumes[1]
    # contains spherical volumes of viperdb viral capsids
    capsid_final = calc_sphericalVol()
    # only ds spherical
    dissected_lists = list_capsidVOL_ds()
    capss_vol = dissected_lists[0]
    new_capsid_final = dissected_lists[1]
    # decide path for later based on comparison between the genome volume size and capsid size
    for g in compare_genomeandcapsidVol():
        if g[1] =='dec1' and g[0] ==True: #only for available capsids, no conditionals for dec2 ? is necessary?
            counter = g[2]
            dec1 = True
            c1_posav = g[3]
        elif g[1] =='dec2' and g[4]== True: #there are some available and also some that are unavailable
            counter = g[2]
            dec2 = True
        elif g[1] =='dec2' and g[4] == False: #there are no available capsid
            counter = g[2]
            dec2 = True
            dec1 = False

    # call in lists created for the capsids
    approvedlists = approvlisting()
    fam_lisF = approvedlists[0]
    fam_lisF2 = approvedlists[1]
    capsVlis = convertpctoVOL()  # contains list of genome volume/pc for approved viral capsids
    # print summary of PC comparison findings
    sortandfindPCvals2()
    addposofFOUNDPC()
    from P2 import pp1_V
    # CONDITIONS FOR FURTHER CALCULATIONS
    find_len1 = len(pp1_V)  # find the number of volumes that can be used
    if find_len1 >= 1:
        havePCatleastone2()
    elif find_len1 == 0:
        noPCone2()
    # tval calculaiton, with available PCs
    if find_len1 >= 1:
        availablePTCVAL2()
    elif find_len1 == 0:
        NOavailablePCTVAL2()
    TVALreadjustmentsall2()

    print("\n\n" + "-" * 180)
    print("\t" * 20 + "||PROTEIN RATIO||")
    print("-" * 180)
    print("\nThe Adeno-associated virus (Parvoviridae Viral Family ) has a triangulation number of 1. It's total number of protein subunits is equalt to 60.")
    print("The capsid is made up of three structural proteins; VP1, VP2, and VP3. These exist in a ratio of: 1 : 1 : 10"
          "\nReadjustments of this capsid to package your gene can be achieved through alteration of this protein ratio.")
    calcRATIO()
    ratioTABLE()

    print("\n\n" + "-" * 180)
    print("\t" * 14 + "||FINDING A RIBOSOMAL BINDING SITE FROM OUR LIST OF CONSTRUCTS||")
    print("-" * 180)

    from P2 import all_liss
    # unpack new ratio lists
    scaldif_full = []
    full_lis = []
    for first in all_liss:
        scaldif_lis = []
        full1 = []
        for second in first:
            secondpos = first.index(second)
            fullvalcapsid = second * 20
            full1.append(fullvalcapsid)
            ogval = og_numfacetot[secondpos]
            scaldif = fullvalcapsid / ogval
            scaldif_lis.append(scaldif)
        scaldif_full.append(scaldif_lis)
        full_lis.append(
            full1)  # list contains the total number fo each subunit protein required for the new capsid (may include more than one list/)
    # calculate the TEs
    all_vp1 = list(zip(*full_lis))[0]
    all_vp2 = list(zip(*full_lis))[1]
    all_vp3 = list(zip(*full_lis))[2]
    all_vp1scal = list(zip(*scaldif_full))[
        0]  # list of scalar values for vp1, compare with preex list for eval to choose
    all_vp2scal = list(zip(*scaldif_full))[1]
    all_vp3scal = list(zip(*scaldif_full))[2]

    # Calculations for Vp1
    listtest2 = []
    for one in all_vp1:  # unpack vp1 values to calculate TEs
        listtest1 = []
        for i in range(5):
            alpha = 1 + Kr1[int(i)] * Px1[int(i)] * 2850 + Kr1[int(i)] * Px1[int(i)] * int(one)
            denom = 2 * Kr1[int(i)] * Px1[int(i)] * int(one)
            num = math.sqrt((alpha ** 2) - 4 * (((Kr1[int(i)]) ** 2)) * (((Px1[int(i)]) ** 2)) * 2850 * int(one))
            numf = alpha - num
            final_te1 = numf / denom
            listtest1.append(final_te1)
        listtest2.append(listtest1)  # list of list of efficiencies for VP1
    # find closest scal value to choose val in vp1 list of lists  (just vp1)
    findposyscal = []
    for scv1 in all_vp1scal:
        selvalscal = min(scal_og, key=lambda x: abs(x - scv1))
        fpscal = scal_og.index(selvalscal)
        findposyscal.append(fpscal)  # contains the list of positions for the selected constructs
    collectselect = []
    collectdescp = []
    for it1, it2 in zip(findposyscal, listtest2):
        findsdescp = lisdescp[int(it1)]
        collectdescp.append(findsdescp)  # contains the list of selected descp for vp1
        selectconstruct = it2[int(it1)]  # find the most appropriate construct within the calc vp1 list
        collectselect.append(selectconstruct)  # contains list of selected constructs for vp1

    # Calculations VP2
    listtest3 = []
    for two in all_vp2:  # unpack vp1 values to calculate TEs
        listtest4 = []
        for i in range(5):
            alpha = 1 + Kr1[int(i)] * Px1[int(i)] * 2850 + Kr1[int(i)] * Px1[int(i)] * int(two)
            denom = 2 * Kr1[int(i)] * Px1[int(i)] * int(one)
            num = math.sqrt((alpha ** 2) - 4 * (((Kr1[int(i)]) ** 2)) * (((Px1[int(i)]) ** 2)) * 2850 * int(two))
            numf = alpha - num
            final_te2 = numf / denom
            listtest4.append(final_te2)
        listtest3.append(listtest4)  # list of list of efficiencies for VP1
    # find closest scal value to choose val in vp1 list of lists  (just vp1)
    findposyscal2 = []
    for scv2 in all_vp2scal:
        selvalscal2 = min(scal_og, key=lambda x: abs(x - scv2))
        fpscal2 = scal_og.index(selvalscal2)
        findposyscal2.append(fpscal2)  # contains the list of positions for the selected constructs
    collectselect2 = []
    collectdescp2 = []
    for it3, it4 in zip(findposyscal2, listtest3):
        findsdescp2 = lisdescp[int(it3)]
        collectdescp2.append(findsdescp2)  # contains the list of selected descp for vp1
        selectconstruct2 = it4[int(it3)]  # find the most appropriate construct within the calc vp1 list
        collectselect2.append(selectconstruct2)  # contains list of selected constructs for vp1

    # Calculation VP3
    listtest5 = []
    for three in all_vp3:  # unpack vp1 values to calculate TEs
        listtest6 = []
        for i in range(5):
            alpha = 1 + Kr1[int(i)] * Px1[int(i)] * 2850 + Kr1[int(i)] * Px1[int(i)] * int(three)
            denom = 2 * Kr1[int(i)] * Px1[int(i)] * int(one)
            num = math.sqrt((alpha ** 2) - 4 * (((Kr1[int(i)]) ** 2)) * (((Px1[int(i)]) ** 2)) * 2850 * int(three))
            numf = alpha - num
            final_te3 = numf / denom
            listtest6.append(final_te3)
        listtest5.append(listtest6)  # list of list of efficiencies for VP1
    # find closest scal value to choose val in vp1 list of lists  (just vp1)
    findposyscal3 = []
    for scv3 in all_vp3scal:
        selvalscal3 = min(scal_og, key=lambda x: abs(x - scv3))
        fpscal3 = scal_og.index(selvalscal3)
        findposyscal3.append(fpscal3)  # contains the list of positions for the selected constructs
    collectselect3 = []
    collectdescp3 = []
    for it5, it6 in zip(findposyscal3, listtest5):
        findsdescp3 = lisdescp[int(it5)]
        collectdescp3.append(findsdescp2)  # contains the list of selected descp for vp1
        selectconstruct3 = it6[int(it5)]  # find the most appropriate construct within the calc vp1 list
        collectselect3.append(selectconstruct3)  # contains list of selected constructs for vp1
    valte1 = ['%.2f' % elem for elem in collectselect]
    valte2 = ['%.2f' % elem for elem in collectselect2]
    valte3 = ['%.2f' % elem for elem in collectselect3]

    new_lst = [list(x) for x in zip(valte1, valte2, valte3)]
    select_top = max(new_lst)  # ratio with the largest total tes
    find_top = new_lst.index(select_top)
    select_topRatio = all_liss[find_top]  # select ratio is list - found

    second_lst = [list(x) for x in zip(collectdescp, collectdescp2,collectdescp3)]
    select_descp = second_lst[find_top]

    print("\n" * 3 + "\t" * 16 + "_Summary of Selected Wetlab Constructs_\n")
    titles = ['New Ratio', 'VP1 Construct', 'TE1', ' VP2 Construct', 'TE2',
              'VP3 Construct', 'TE3']
    data = [titles] + list(zip(all_liss, collectdescp, valte1, collectdescp2, valte2, collectdescp3, valte3))
    for i, d in enumerate(data):
        line = ''.join(str(x).ljust(23) for x in d)
        print("\t" + line)
        if i == 0:
            print('-' * len(line) + '\n')

    print("\n\nDescription of Outputs")
    print("-" * 22)
    print(
        "TE - Translation Efficiency\n*Calculations for model uses mRNA-folding dynamics and ribosome-binding dynamics,\nPercentage values (ie. x100) reflect the probability of the ribosome binding to VP1, VP2, or VP3 mRNA transcripts.")
    print("\nP1 - Strong Promoter\nP2 - Medium Promoter\nRBS(1/2/3) - Weak, Medium, Strong - Ribosomal Binding Sites")

    print("\n\n" + "-" * 180)
    print("\t" * 22 + "||END||")
    print("-" * 180)
    print('\n' * 2 + '\t' * 12 + '--------- Final Selection of Suitable Wet Construct and Protein Ratio ----------\n')
    print('\t' * 12 + "Selected Protein Ratio: " + ' : '.join(map(str, select_topRatio)))
    print('\n' + '\t' * 12 + "Translation Efficiencies of VP1, VP2, VP3 Subunit Proteins: " + ', '.join(
        map(str, select_top)))
    print('\n' + '\t' * 12 + "Wet Construct Selection: " + ', '.join(map(str, select_descp)))
    print('\n' + '\t' * 12 + '-' * 85)

    #lauch page for graph

    x = select_topRatio
    y = select_top
    z = select_descp
    prot = ['VP1', 'VP2', 'VP3']
    print(select_top)
    print(prot)

    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        row_heights=[0.4, 0.8],
        specs=[[{"type": "table"}],
               [{"type": "scatter"}]]
    )

    fig.add_trace(go.Bar(x=['VP1', 'VP2', 'VP3'],
                         y=y),
                  row=2, col=1
                  )

    fig.add_trace(
        go.Table(
            header=dict(
                values=["AAV Capsid Subunit Protein", "Selected Protein Ratio Value", "Translation Efficiency", "Selected Wet Lab Construct"],
                font=dict(size=13),
                align="center"
            ),
            cells=dict(
                values=[prot, x, y, z],
                align="center")
        ),
        row=1, col=1
    )
    fig.update_layout(
        height=800,
        showlegend=False,
        title_text="Translation of Efficiency of Selected Wet Lab Construct",

    )

    fig.show()

    # this restart the script again for another gene input -----end with this!
    restart = input("\n\n\nWould you like enter another gene length?")
    if restart == "yes" or restart == "y":
        script()
    elif restart == "n" or restart == "no":
        print("Goodbye.")


script()
