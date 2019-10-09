'''
CapsidOptimizer

Calculate new icosahedral geometries to design novel theoretical viral capsids with optimised 
packaging capacities for the delivery of a select gene associated with a rare genetic disease. 

'''
from P2 import*
def script():
    table_presets()
    gene = ask_user()
    if gene not in gen_sizeR:
        print("There is no gene length that matches this in the table, please try again.")
        script()
    if gene == 0:  # if the user inputs 0, should restart script
        print("Your gene length must be greater than 0. Please try again.\n")
        script()

    #ds and ss genome volumes - cylindrical volumes
    gene_volumes = gene_vol()
    gen_vol = gene_volumes[0]
    gen_volss = gene_volumes[1]
    #contains spherical volumes of viperdb viral capsids
    capsid_final = calc_sphericalVol()
    #only ds spherical
    dissected_lists = list_capsidVOL_ds()
    capss_vol = dissected_lists[0]
    new_capsid_final = dissected_lists[1]
    #decide path for later based on comparison between the genome volume size and capsid size
    for g in compare_genomeandcapsidVol():
        if g[1] =='dec1' and g[0] ==True: #only for available capsids, no conditionals for dec2 ? is necessary?
            counter = g[2]
            dec1 = True
        elif g[1] =='dec2' and g[4]== True: #there are some available and also some that are unavailable
            counter = g[2]
            dec2 = True
        elif g[1] =='dec2' and g[4] == False: #there are no available capsid
            counter = g[2]
            dec2 = True
            dec1 = False
    if counter == 0:
        nocapsidsVol()
    if dec1 == True and dec2 == False:
        allcapsidsVOL()
    if dec1 == True and dec2 == True:
        somecapsidsVOL()


    # call in lists created for the capsids
    approvedlists = approvlisting()
    fam_lisF = approvedlists[0]
    fam_lisF2 = approvedlists[1]

    # Start Section for Determination of Capsid Volumes
    print("-" * 180)
    print("\t" * 18 + "||CAPSID VOLUMES OF VIRAL VECTORS||")
    print("-" * 180)
    # print the dimension genome volumes - cylindrical dimensions
    print("\nDouble Stranded Genome Volume (A): " + str(gen_vol))
    print("\nSingle Stranded Genome Volume (A): " + str(gen_volss))

    print("\n\nWHAT ARE THE VIRUSES CURRENTLY BEING USED OR TESTED FOR THERAPEUTICS?")
    print("\nViral families currently being researched for viral delivery: " + ", ".join(fam_lisF))
    print("\nViral families that have been approved for the delivery of therapeutics: " + ", ".join(fam_lisF2))

    #Find the approved viral vectors
    cond_2 = findapprovedcapsids()
    findaprovTABLE()
    if cond_2 == True:
        someavailableSUM()
    if cond_2 ==False:
        noavailableSUM()

    #PACKAGING CAPACITITES
    print("-" * 180)
    print("\t" * 18 + "||PACKAGING CAPACITY OF VIRAL VECTORS||")
    print("-" * 180)
    print("*evaluation of packing capacity is achieved by a comparison between the 'actual' nucleic acid volumes between the selected gene and allowable capacity of each viral capsid*")
    print("-this assessment also considers the type of genome each viral capsids perferentially packages (ss versus ds nucleic acid)")

    capsVlis = convertpctoVOL() #contains list of genome volume/pc for approved viral capsids

    #print summary of PC comparison findings
    sortandfindPCvals()
    addposofFOUNDPC()
    from P2 import pp1_V
    # CONDITIONS FOR FURTHER CALCULATIONS
    find_len1 = len(pp1_V)  # find the number of volumes that can be used
    if find_len1 >= 1:
        havePCatleastone()
    elif find_len1 == 0:
        noPCone()

    #Determine the new Triangulation Number
    print("\n\n" + "-" * 180)
    print("\t" * 18 + "||TRIANGULATION NUMBER||")
    print("-" * 180)

    #tval calculaiton, with available PCs
    decide_rbs1 = False
    decide_rbs2 = False
    if find_len1 >= 1:
        availablePCTVAL()
        from P2 import decide_rbs1
    elif find_len1 == 0:
        NOavailablePCTVAL()
        from P2 import decide_rbs2
    #create final summary table of total Tval adjustments
    TVALreadjustmentsall()


    #Determine the Capsid Stability - else also add additional biophysic calculations?
    print("\n\n" + "-" * 180)
    print("\t" * 18 + "||EVALUATION OF CAPSID STABILITY||")
    print("-" * 180)
    #verdict dependent values - capsid stab
    if decide_rbs1 == True:
        capsidstabwithAVAL()
    if decide_rbs2 == True:
        capsidstabwithNOOAVAL()
    '''
    add calculations for biophysics  - include final summary table for final capsid selection 
    '''
    #final table for current therapeutic using viral capsids
    therapeuticsTABLE()

    print("\n\n" + "-" * 180)
    print("\t" * 22 + "||END||")
    print("-" * 180)

    # print summary tables of final output
    if find_len1 >= 1:
        FINALTABLE1()
    elif find_len1 == 0:
        FINALTABLE2()



    if find_len1 >=1:
        conv_vol_diamAVCAPS()
        summaryplotall()
    elif find_len1 ==0:
        conv_vol_diamNOCAPS()
        summaryplotall()



    #this restart the script again for another gene input -----end with this!
    restart = input("\n\n\nWould you like enter another gene length?")
    if restart == "yes" or restart == "y":
        script()
    elif restart == "n" or restart == "no":
        print("Goodbye.")
script()



















