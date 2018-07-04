from collections import defaultdict
import Dedup

def get_prediction(eq_obj):
    predictions = {}
    for gene in eq_obj.genes_list:
        predictions[gene] = 0

    # extract gene-level umis
    gene_umis = defaultdict(lambda : defaultdict(int))
    for eqId, label in enumerate(eq_obj.labels):
        gene_set = set([])
        for txp in label:
            txp_name = eq_obj.txps_list[txp]
            gene = eq_obj.txp_to_gene_dict[txp_name]
            gene_set.add(gene)
        if len(gene_set) == 1:
            gene_name = list(gene_set)[0]
            for umi, count in eq_obj.umis[eqId].items():
                gene_umis[gene_name][umi] += count
        else:
            print ("ERROR: gene set for eqclass {} not 1 length for Utools".format(label))

    print ("Utools")
    print (gene_umis)
    for gene, umis in gene_umis.items():
        predictions[gene] = Dedup.dedup(umis, True)

    return predictions
