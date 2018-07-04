def read_quant_alevin(base):
    alv = pd.read_table(base+"quants_mat.csv", sep=",", header=None)
    names = pd.read_table(base+"quants_mat_rows.txt", header=None)
    genes = pd.read_table(base+"quants_mat_cols.txt", header=None)
    alv.drop([len(alv.columns)-1], axis=1, inplace=True)
    alv.columns = genes[0].values
    alv.index = names[0].values
    alv = alv.loc[:, (alv != 0).any(axis=0)]
    return alv

def read_tenx(genome, base):
    matrices_dir = base
    human_matrix_dir = os.path.join(matrices_dir, genome)
    mat = scipy.io.mmread(os.path.join(human_matrix_dir, "matrix.mtx")).toarray()

    genes_path = os.path.join(human_matrix_dir, "genes.tsv")
    gene_ids = [row[0] for row in csv.reader(open(genes_path), delimiter="\t")]
    gene_names = [row[1] for row in csv.reader(open(genes_path), delimiter="\t")]

    barcodes_path = os.path.join(human_matrix_dir, "barcodes.tsv")
    barcodes = [row[0][:-2] for row in csv.reader(open(barcodes_path), delimiter="\t")]

    cr = pd.DataFrame(mat).T
    cr.index = barcodes
    cr.columns = gene_ids

    return cr

def read_naive(infile):
    #utools import
    naive = pd.read_table(infile, index_col=0)
    return naive

def read_sf_bin(bf):
    gene_names = pd.read_table(bf+"/quants_mat_cols.txt", header=None)
    cell_names = pd.read_table(bf+"/quants_mat_rows.txt", header=None)

    numTxps = len(gene_names)
    header_struct = struct.Struct("d"*numTxps)
    with gzip.open(bf) as f:
        count = 0
        no_read_count = 0
        tot_read_count = 0
        umiCounts = []
        while True:
            count += 1
            if count%100 == 0:
                print ("\r Done reading " + str(count) + " cells", end="")
                sys.stdout.flush()
            try:
                cell_counts = header_struct.unpack_from( f.read(header_struct.size) )
            except:
                print ("\nRead total " + str(count-1) + " cells")
                print ("Found total " + str(tot_read_count) + " reads")
                break
            read_count = 0.0
            for x in cell_counts:
                read_count += float(x)
            tot_read_count += read_count
            if read_count > 0.0:
                umiCounts.append( cell_counts )
            else:
                no_read_count += 1
    print ("No Read Count Cells: "+str(no_read_count))
    df = pd.DataFrame(umiCounts)
    df.index = cell_names
    df.columns = gene_names

    return df
