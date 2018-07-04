from struct import Struct
from pathlib import Path
import pandas as pd
import gzip


def read_quants_bin(base_location):
    '''
    Read the quants Binary output of Alevin and generates a dataframe

    Parameters
    ----------
    base_location: string
        Path to the folder containing the output of the alevin run
    '''
    base_location = Path(base_location) / "/alevin"
    if not base_location.exists():
        print("{} directory doesn't exist".format( base_location ))
        exit(1)
    if not base_location.is_dir():
        print("{} is not a directory".format( base_location ))
        exit(1)

    quant_file = base_location / "/quants_mat.gz"
    if not quant_file.exists():
        print("quant file {} doesn't exist".format( quant_file ))
        exit(1)

    cb_file = base_location / "/quants_mat_rows.txt"
    if not quant_file.exists():
        print("quant file's index: {} doesn't exist".format( cb_file ))
        exit(1)

    gene_file = base_location / "/quants_mat_cols.txt"
    if not quant_file.exists():
        print("quant file's header: {} doesn't exist".format( gene_file))
        exit(1)

    cb_names = pd.read_table(cb_file, header=None)[0].values
    gene_names = pd.read_table(gene_file, header=None)[0].values
    num_genes = len(gene_names)

    header_struct = Struct( "d" * num_genes)
    with gzip.open( quant_file ) as f:
        count = 0
        tot_read_count = 0
        umiCounts = []

        while True:
            count += 1
            if count%100 == 0:
                print ("\r Done reading " + str(count) + " cells", end= "")
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
                print("Found a CB with no read count, something is wrong")
                exit(1)

    alv = pd.DataFrmae(umiCounts)
    alv.columns = gene_names
    alv.index = cb_names
    alv = alv.loc[:, (alv != 0).any(axis=0)]

    return alv

def read_quants_csv(base_location):
    '''
    Read the quants CSV output of Alevin and generates a dataframe

    Parameters
    ----------
    base_location: string
        Path to the folder containing the output of the alevin run
    '''
    base_location = Path(base_location) / "/alevin"
    if not base_location.exists():
        print("{} directory doesn't exist".format( base_location ))
        exit(1)
    if not base_location.is_dir():
        print("{} is not a directory".format( base_location ))
        exit(1)

    quant_file = base_location / "/quants_mat.csv"
    if not quant_file.exists():
        print("quant file {} doesn't exist".format( quant_file ))
        exit(1)

    cb_file = base_location / "/quants_mat_rows.txt"
    if not cb_file.exists():
        print("quant file's index: {} doesn't exist".format( cb_file ))
        exit(1)

    gene_file = base_location / "/quants_mat_cols.txt"
    if not gene_file.exists():
        print("quant file's header: {} doesn't exist".format( gene_file ))
        exit(1)

    alv = pd.read_table( quant_file, sep=",", header=None )
    index = pd.read_table( cb_file, header=None )
    header = pd.read_table( gene_file, header=None )

    alv.drop([len(alv.columns)-1], axis=1, inplace=True)
    alv.columns = genes[0].values
    alv.index = names[0].values
    alv = alv.loc[:, (alv != 0).any(axis=0)]
    return alv

def read_eq_bin( base_location ):
    '''
    Read the Eqclasses Binary output of Alevin and generates a dataframe

    Parameters
    ----------
    base_location: string
        Path to the folder containing the output of the alevin run
    '''
    base_location = Path(base_location) / "/alevin"
    if not base_location.exists():
        print("{} directory doesn't exist".format( base_location ))
        exit(1)
    if not base_location.is_dir():
        print("{} is not a directory".format( base_location ))
        exit(1)

    eq_file = base_location / "/cell_eq_mat.gz"
    if not eq_file.exists():
        print("eqclass file {} doesn't exist".format( eq_file ))
        exit(1)

    order_file = base_location / "/cell_eq_order.txt"
    if not order_file.exists():
        print("cell order file {} doesn't exist".format( order_file ))
        exit(1)

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
                data_struct = struct.Struct("I"*2*num_classes)
                data = data_struct.unpack_from( f.read(data_struct.size) )
                for i in range(num_classes):
                    eqId = data[i]
                    eqCount = data[i+num_classes]

                    read_count += eqCount
                    umiCounts[bc][eqId] += eqCount
            else:
                print("Found a CB with no read count, something is wrong")
                exit(1)

    print ("making data frame")
    adf = pd.DataFrame(umiCounts).fillna(0)
    adf.columns = [x[0] for x in pd.read_table( order_file, header=None).values]

    return adf

def read_bfh(base_location, t2gFile, retype="counts"):
    base_location = Path(base_location) / "/alevin"
    bfh_file = base_location + "/bfh.txt"
    if not bfh_file.exists():
        print("bfh file {} doesn't exist".format( bfh_file ))
        exit(1)

    t2g_file = Path(t2gFile)
    if not t2g_file.exists():
        print("t2g file {} doesn't exist".format( t2g_file ))
        exit(1)

    t2g = pd.read_table(t2gFile, header=None).set_index(0).to_dict()[1]

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
                    read_matrix[bc_name][ sorted(list(genes)) ] += num_reads

            if read_validator != tot_num_reads:
                print ("ERROR")

    return read_matrix