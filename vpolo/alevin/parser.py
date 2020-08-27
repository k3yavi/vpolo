from __future__ import print_function
from collections import defaultdict
from struct import Struct
import numpy as np
import pandas as pd
import gzip
import sys
import os
import sce
from scipy.io import mmread
from scipy.sparse import csr_matrix
def read_quants_bin(base_location, clipped=False, mtype="data", rmode="rust"):
    '''
    Read the quants Sparse Binary output of Alevin and generates a dataframe
    Parameters
    ----------
    base_location: string
        Path to the folder containing the output of the alevin run
    clipped: bool (default False)
        Clip off all zero rows and columns
    mtype: "[data(default), tier, var, mean]"
        Alevin's matrix type to load into memory
    '''
    if not os.path.isdir(base_location):
        print("{} is not a directory".format( base_location ))
        sys.exit(1)

    base_location = os.path.join(base_location, "alevin")
    
    print(base_location)
    if not os.path.exists(base_location):
        print("{} directory doesn't exist".format( base_location ))
        sys.exit(1)

    data_type = "f"
    if mtype == "data":
        quant_file = os.path.join(base_location, "quants_mat.gz")
    elif mtype == "tier":
        data_type = "B"
        quant_file = os.path.join(base_location, "quants_tier_mat.gz")
    elif mtype == "mean":
        quant_file = os.path.join(base_location, "quants_mean_mat.gz")
    elif mtype == "var":
        quant_file = os.path.join(base_location, "quants_var_mat.gz")
    else:
        print("wrong mtype:".format( mtype ))
        sys.exit(1)

    if not os.path.exists(quant_file):
        print("quant file {} doesn't exist".format( quant_file ))
        sys.exit(1)

    if mtype in ["mean", "var"]:
        cb_file = os.path.join(base_location, "quants_boot_rows.txt")
    else:
        cb_file = os.path.join(base_location, "quants_mat_rows.txt")

    if not os.path.exists(cb_file):
        print("quant file's index: {} doesn't exist".format( cb_file ))
        sys.exit(1)

    gene_file = os.path.join(base_location, "quants_mat_cols.txt")
    
    if not os.path.exists(gene_file):
        print("quant file's header: {} doesn't exist".format(gene_file))
        sys.exit(1)

    cb_names = pd.read_csv(cb_file, header=None)[0].values
    gene_names = pd.read_csv(gene_file, header=None)[0].values
    num_genes = len(gene_names)
    num_cbs = len(cb_names)
    num_entries = int(np.ceil(num_genes/8))

    if rmode == "rust":
        print("Using rust mode with {} rows and {} columns".format(num_cbs, num_genes))
        mat = sce.read_quants(quant_file, num_cbs, num_genes)
        umi_matrix = csr_matrix((mat[2], mat[1], mat[0]), shape=(num_cbs, num_genes)).todense()
    else:
        with gzip.open( quant_file ) as f:
            line_count = 0
            tot_umi_count = 0
            umi_matrix = []

            header_struct = Struct( "B" * num_entries)
            while True:
                line_count += 1
                if line_count%100 == 0:
                    print ("\r Done reading " + str(line_count) + " cells", end= "")
                    sys.stdout.flush()
                try:
                    num_exp_genes = 0
                    exp_counts = header_struct.unpack_from( f.read(header_struct.size) )
                    for exp_count in exp_counts:
                        num_exp_genes += bin(exp_count).count("1")

                    data_struct = Struct( data_type * num_exp_genes)
                    sparse_cell_counts_vec = list(data_struct.unpack_from( f.read(data_struct.size) ))[::-1]
                    cell_umi_counts = sum(sparse_cell_counts_vec)

                except:
                    print ("\nRead total " + str(line_count-1) + " cells")
                    print ("Found total " + str(tot_umi_count) + " reads")
                    break

                if cell_umi_counts > 0.0:
                    tot_umi_count += cell_umi_counts

                    cell_counts_vec = []
                    for exp_count in exp_counts:
                        for bit in format(exp_count, '08b'):
                            if len(cell_counts_vec) >= num_genes:
                                break

                            if bit == '0':
                                cell_counts_vec.append(0.0)
                            else:
                                abund = sparse_cell_counts_vec.pop()
                                cell_counts_vec.append(abund)

                    if len(sparse_cell_counts_vec) > 0:
                        print("Failure in consumption of data")
                        print("left with {} entry(ies)".format(len(sparse_cell_counts_vec)))
                    umi_matrix.append( cell_counts_vec )
                else:
                    print("Found a CB with no read count, something is wrong")
                    sys.exit(1)


    alv = pd.DataFrame(umi_matrix)
    alv.columns = gene_names
    alv.index = cb_names
    if clipped:
        alv = alv.loc[:, (alv != 0).any(axis=0)]

    return alv

def read_fry_bootstraps_bin(
    base_location,
    num_bootstraps,
):
    '''
    Read the bootstraps as EDS output of Alevin and generates a dictionary of 
    dataframes, one for each cell. The dimension of each dataframe is 
    num_bootstraps x num_genes
    Parameters
    ----------
    base_location: string
        Path to the folder containing the output of the alevin run
    num_bootstraps: int (default False)
        number of bootstraps
    density: "[sparse (default), dense ]" (not in use)
        Load sparse alevin output or dense output(<v0.14.0)
    mtype: "[data(default), tier, var, mean]" (not in use)
        Alevin's matrix type to load into memory

    '''
    if not os.path.isdir(base_location):
        print("{} is not a directory".format( base_location ))
        sys.exit(1)

    #base_location = os.path.join(base_location, "alevin")
    print(base_location)
    if not os.path.exists(base_location):
        print("{} directory doesn't exist".format( base_location ))
        sys.exit(1)

    data_type = "f"
    quant_file = os.path.join(base_location, "bootstraps.eds.gz")

    if not os.path.exists(quant_file):
        print("bootstrap file {} doesn't exist".format( quant_file ))
        sys.exit(1)
        
    cb_file = os.path.join(base_location, "barcodes.txt")
    if not os.path.exists(cb_file):
        print("quant file's index: {} doesn't exist".format( cb_file ))
        sys.exit(1)
    gene_file = os.path.join(base_location, "gene_names.txt")
    if not os.path.exists(gene_file):
        print("quant file's header: {} doesn't exist".format( gene_file))
        sys.exit(1)

    cb_names = pd.read_csv(cb_file, header=None, sep='\t')[1].values
    gene_names = pd.read_csv(gene_file, header=None)[0].values
    num_genes = len(gene_names)
    num_entries = int(np.ceil(num_genes/8))
    
    bootstrap_dict = {}
    cb_index = 0
    with gzip.open( quant_file ) as f:
        line_count = 0
        tot_umi_count = 0
        umi_matrix = []

        header_struct = Struct( "B" * num_entries)

        while True:
            line_count += 1
            if line_count%100 == 0:
                print ("\r Done reading " + str(line_count) + " cells", end= "")
                sys.stdout.flush()
                
            if line_count%num_bootstraps == 0:
                alv = pd.DataFrame(umi_matrix)
                alv.columns = gene_names
                bootstrap_dict[cb_names[cb_index]] = alv
                cb_index += 1
                umi_matrix = []
                
            try:
                num_exp_genes = 0
                exp_counts = header_struct.unpack_from( f.read(header_struct.size) )
                for exp_count in exp_counts:
                    num_exp_genes += bin(exp_count).count("1")

                data_struct = Struct( data_type * num_exp_genes)
                sparse_cell_counts_vec = list(data_struct.unpack_from( f.read(data_struct.size) ))[::-1]
                cell_umi_counts = sum(sparse_cell_counts_vec)

            except:
                print ("\nRead total " + str(line_count-1) + " cells")
                print ("Found total " + str(tot_umi_count) + " reads")
                break

            if cell_umi_counts > 0.0:
                tot_umi_count += cell_umi_counts

                cell_counts_vec = []
                for exp_count in exp_counts:
                    for bit in format(exp_count, '08b'):
                        if len(cell_counts_vec) >= num_genes:
                            break

                        if bit == '0':
                            cell_counts_vec.append(0.0)
                        else:
                            abund = sparse_cell_counts_vec.pop()
                            cell_counts_vec.append(abund)

                if len(sparse_cell_counts_vec) > 0:
                    print("Failure in consumption of data")
                    print("left with {} entry(ies)".format(len(sparse_cell_counts_vec)))
                umi_matrix.append( cell_counts_vec )
            else:
                print("Found a CB with no read count, something is wrong")
                sys.exit(1)
        
    #alv.index = cb_names
    return bootstrap_dict


def read_eq_bin( base_location ):
    '''
    Read the Eqclasses Binary output of Alevin and generates a dataframe

    Parameters
    ----------
    base_location: string
        Path to the folder containing the output of the alevin run
    '''
    base_location = os.path.join(base_location, "alevin")
    if not os.path.isdir(base_location):
        print("{} is not a directory".format( base_location ))
        sys.exit(1)

    eq_file = os.path.join(base_location, "cell_eq_mat.gz")
    if not os.path.exists(eq_file):
        print("eqclass file {} doesn't exist".format( eq_file ))
        sys.exit(1)

    order_file = os.path.join(base_location, "cell_eq_order.txt")
    if not os.path.exists(order_file):
        print("cell order file {} doesn't exist".format( order_file ))
        sys.exit(1)

    header_struct = Struct("Q"*2)
    with gzip.open( eq_file ) as f:
        count = 0
        read_count = 0
        umiCounts = defaultdict(lambda: defaultdict(int))
        while True:
            count += 1
            if count%100 == 0:
                print ("\r Done reading " + str(count) + " cells", end= "")
                sys.stdout.flush()
            try:
                bc, num_classes = header_struct.unpack_from( f.read(header_struct.size) )
            except:
                print ("\nRead total " + str(count-1) + " cells")
                print ("Found total " + str(read_count) + " reads")
                break

            if num_classes != 0:
                data_struct = Struct("I"*2*num_classes)
                data = data_struct.unpack_from( f.read(data_struct.size) )
                for i in range(num_classes):
                    eqId = data[i]
                    eqCount = data[i+num_classes]

                    read_count += eqCount
                    umiCounts[bc][eqId] += eqCount
            else:
                print("Found a CB with no read count, something is wrong")
                sys.exit(1)

    print ("making data frame")
    adf = pd.DataFrame(umiCounts).fillna(0)
    adf.columns = [x[0] for x in pd.read_csv( order_file, header=None).values]

    return adf

def read_arborescences(base_location):
    '''
    Read the quants Sparse Binary output of Alevin and generates a dataframe
    Parameters
    ----------
    base_location: string
        Path to the folder containing the output of the alevin run
    Returns:
        a mult-level dictionary with CB -> Gene -> Arborescence_length -> Frequency
    '''
    if not os.path.isdir(base_location):
        print("{} is not a directory".format( base_location ))
        sys.exit(1)

    base_location = os.path.join(base_location, "alevin")
    print(base_location)
    if not os.path.exists(base_location):
        print("{} directory doesn't exist".format( base_location ))
        sys.exit(1)

    arbo_file = os.path.join(base_location, "arborescence_dump.txt.gz")
    if not os.path.exists(arbo_file):
        print("quant file {} doesn't exist".format( quant_file ))
        sys.exit(1)

    cell_gene_arbo_data = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    with gzip.open( arbo_file ) as f:
        cell_count = 0
        total_fragments = 0

        while True:
            cell_count += 1
            if cell_count % 10 == 0:
                print ("\r Done reading " + str(cell_count) + " cells", end= "")
                sys.stdout.flush()
                
            # each cell has a header line of CB name, number of expressed genes, 
            # total number of fragments
            try:
                (cell_name, num_exp_genes, num_fragments) = f.readline().decode().strip().split()
                total_fragments += int(num_fragments)
            except:
                print("\nRead Total {} cells w/ {} fragments".format(cell_count-1, total_fragments))
                break
            
            # each Cell is followed by gene with at least one mapped fragment
            # and has gene_id, nuumber_of_types_of arbos, 
            # [#frags in arbo, frequency of such arbo]+
            for _ in range(int(num_exp_genes)):
                (gid, num_arbo_types, *arbo_data) = f.readline().decode().strip().split()
                gid = int(gid)
                
                # iterating over each gene's arborescence data
                for i in range(int(num_arbo_types)):
                    arbo_index = 2*i
                    arbo_length = int(arbo_data[arbo_index])
                    arbo_frequency = int(arbo_data[arbo_index+1])
                    
                    cell_gene_arbo_data[cell_name][gid][arbo_length] += arbo_frequency

    return cell_gene_arbo_data

def read_bfh(base_location, t2g_file, retype="counts"):
    base_location = os.path.join(base_location, "alevin")
    if not os.path.isdir(base_location):
        print("{} is not a directory".format( base_location ))
        sys.exit(1)

    bfh_file = os.path.join(base_location, "bfh.txt")
    if not os.path.exists(bfh_file):
        print("bfh file {} doesn't exist".format( bfh_file ))
        sys.exit(1)

    if not os.path.exists(t2g_file):
        print("t2g file {} doesn't exist".format( t2g_file ))
        sys.exit(1)

    t2g = pd.read_csv(t2g_file, header=None, sep="\t").set_index(0).to_dict()[1]

    if retype == "counts":
        read_matrix = defaultdict(lambda : defaultdict(int))
    elif retype == "umis":
        read_matrix = defaultdict(lambda : defaultdict(lambda : defaultdict(int)))

    # read in Alevin
    with open( bfh_file ) as f:
        T = int(f.readline())
        C = int(f.readline())
        E = int(f.readline())

        tname = []
        bc_id_to_name = []

        for _ in range(T):
            tname.append(f.readline().strip())

        for _ in range(C):
            bc_id_to_name.append(f.readline().strip())

        for idx,line in enumerate(f):
            toks = line.strip().split()
            num_labels = int(toks[0])
            tot_num_reads = int(toks[num_labels+1])
            genes = set([])

            for txp in toks[1:num_labels+1]:
                genes.add(t2g[ tname[int(txp)] ])

            idx = num_labels+2
            num_bcs = int(toks[idx])
            read_validator = 0

            for _ in range(num_bcs):
                idx += 1
                bc_name = bc_id_to_name[ int(toks[idx]) ]
                idx += 1
                num_umi = int(toks[idx])
                num_reads = 0

                for _ in range(num_umi):
                    idx += 2
                    num_reads += int(toks[idx])
                    if retype == "umis":
                        read_matrix[bc_name][ tuple(sorted(toks[1:num_labels+1])) ][ toks[idx-1] ] += int(toks[idx])

                read_validator += num_reads
                if retype == "counts":
                    read_matrix[bc_name][ tuple(sorted(list(genes))) ] += num_reads

            if read_validator != tot_num_reads:
                print ("ERROR")

    return read_matrix

def read_tenx(base, version=2):
    '''
    Specify the path to the folder containing matrix.mtx file
    '''
    if version == 2:
        mat = mmread(os.path.join(base, "matrix.mtx")).toarray()

        genes_path = os.path.join(base, "genes.tsv")
        genes = pd.read_csv(genes_path, header=None)[0].values

        barcodes_path = os.path.join(base, "barcodes.tsv")
        barcodes = pd.read_csv(barcodes_path, header=None)[0].values
    elif version == 3:
        mat_file = os.path.join(base, "matrix.mtx.gz")
        with gzip.open(mat_file) as f:
            mat = mmread(f).toarray()

        genes_path = os.path.join(base, "features.tsv.gz")
        with gzip.open(genes_path) as f:
            genes = pd.read_csv(f, header=None)[0].values

        barcodes_path = os.path.join(base, "barcodes.tsv.gz")
        with gzip.open(barcodes_path) as f:
            barcodes = pd.read_csv(f, header=None)[0].values
    else:
        print("Wrong version")

    cr = pd.DataFrame(mat).T
    cr.index = [x.strip().split("-")[0] for x in barcodes]
    cr.columns = genes

    return cr

def read_umi_tools(infile):
    '''
    Specify the umi_tools count output file
    '''
    naive = pd.read_csv(infile, index_col=0, sep="\t")
    return naive.T

def read_umi_graph(base_location, out_location, kind="dot"):
    '''
    A function to read the per cell level UMI graph output from Alevin
    i.e. a file with name cel_umi_graphs.gz and dumps per cell level
    separate dot(viz) file
    '''
    if not os.path.isdir(base_location):
        print("{} directory doesn't exist".format(base_location))
        sys.exit(1)

    if not os.path.isdir(out_location):
        print("{} directory doesn't exist".format(out_location))
        sys.exit(1)

    base_location = os.path.join(base_location, "alevin")
    if not os.path.isdir(base_location):
        print("{} is not a directory".format(base_location))
        sys.exit(1)

    graph_file = os.path.join(base_location, "cell_umi_graphs.gz")
    if not os.path.exists(graph_file):
        print("graph file {} doesn't exist".format(graph_file))
        sys.exit(1)

    with gzip.open(graph_file) as file_handle:
        for cell_index,cell_graph in enumerate(file_handle):
            toks = cell_graph.decode().strip().split("\t")
            fname = toks[0].strip()
            if kind == "dot":
                write_dot(toks[1:], os.path.join(out_location, fname+".dot.gz"))
            print("\r processed {} cell".format(cell_index), end="")

def write_dot(toks, file_name):
    '''
    write per cell level dot file for each cell separately
    '''
    with gzip.open(file_name, 'wb') as f:
        cell_graph = "digraph {} {{".format(file_name)

        # adding vertices
        for vid in range(int(toks[0])):
            cell_graph += "\n" + str(vid)

        #adding edges
        for edge_group in toks[1:]:
            edges = edge_group.strip().split(",")
            cell_graph += "\n{} -> {{ ".format(edges[0])
            cell_graph += ' '.join(edges[1:])
            cell_graph += " }"

        cell_graph += "\n}"

        f.write(cell_graph.encode())
