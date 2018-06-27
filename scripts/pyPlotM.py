import matplotlib.pyplot as plt

def makePlots(tsvFileN, pdfFileN, show):
    #Reading file
    with open(tsvFileN, "r") as ins:
        lines = []
        for line in ins:
            lines.append(line.split('\t'))

    #Calculating maximum number of plots
    maxPlots = 0
    for x in range(1, len(lines)):
        idc = int(lines[x][0]) + 1
        if maxPlots < idc:
            maxPlots = idc

    # Set up the matplotlib figure
    cols = 3
    rows = maxPlots/3
    plt.subplots(3, rows, figsize=(8, 6), sharex=True)

    #Creating array to store the values
    array = []

    for i in range(0, maxPlots):
        arr2 = []
        for j in range(0, 5):
            arr2.append(0)
        array.append(arr2)

    #Filling the array with the values
    for x in range(1, len(lines)):
        idc = int(lines[x][0])
        cnt = int(lines[x][3])
        typ =     lines[x][4].strip()
        nam =     lines[x][1].strip()

        pos = 3;
        if typ == 'synthetic':
            pos = 2
        if typ == 'bcbio':
            pos = 0
        if typ == 'mirge':
            pos = 1
	
        array[idc][pos] = cnt
        array[idc][4]   = nam

    #Plotting the graphs
    plt.figure(1)
    #plt.xlabel('tool')
    #plt.ylabel('Counts')

    p = []
    p.append(array[0][0])
    p.append(array[0][1])
    p.append(array[0][2])
    n = array[0][4]

    for i in range(0, maxPlots):
        del(p[2])
        del(p[1])
        del(p[0])
        p.append(array[i][0])
        p.append(array[i][1])
        p.append(array[i][2])
        n = array[i][4]

        pcd = rows * 100 + cols * 10 + 1 + i
        plt.subplot(pcd)
    
        ax = plt.gca()
        ax.set_facecolor('lightgray')
    
        plt.xticks([1,2,3], ('bcbio', 'mirge', 'synthetic'))
        plt.yticks([0,10,20,30,40,50])
        plt.tick_params(axis='both', which='major', labelsize=8)
        plt.tick_params(axis='both', which='minor', labelsize=8)
    
        plt.bar([1,2,3], p, color='gray')
        plt.title(n)
    
        for i, v in enumerate(p):
            plt.text(i+0.9, 0, str(v), color='black', fontsize='8', fontweight='bold')

    plt.subplots_adjust(top=0.92, bottom=0.10, left=0.10, right=0.95, hspace=0.50, wspace=0.35)
    plt.savefig(pdfFileN, format="pdf")
    if show == 1:
        plt.show()


makePlots("../data/examples/plot/example_count.tsv", "kk.pdf", 1)
