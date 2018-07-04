def deduplicate(eq_obj, predictions):
    for eqId, label in enumerate(eq_obj.labels):
        gene_set = set([])
        for txp in label:
            txp_name = eq_obj.txps_list[txp]
            gene_name = eq_obj.txp_to_gene_dict[txp_name]

            # Add the gene to list of genes found for this eqclass
            gene_set.add(gene_name)
        if len(gene_set) == 1:
            umi_count = len(set(eq_obj.umis[eqId]))
            predictions[list(gene_set)[0]] += umi_count
        else:
            print ("ERROR: Eqclass {} has more than 1 gene and can bias kallisto deduplication"\
                   .format(eqId+1))
            exit(1)

def get_prediction(eq_obj):
    predictions = {}
    for gene in eq_obj.genes_list:
        predictions[gene] = 0

    deduplicate(eq_obj, predictions)

    return predictions