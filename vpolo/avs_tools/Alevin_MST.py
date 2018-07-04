from collections import defaultdict
from itertools import combinations
import Dedup


def get_max_span_networks(gene_net):
    nodes = gene_net.nodes
    edges = {k: v for k, v in gene_net.edges.items()}
    num_networks = len(nodes)

    # breaking ties randomly and very bad implementation
    # basically resolve transitive cases
    for child, parent in edges.items():
        node = list(parent)[0]
        if len(parent) == 1 and node in edges and len(edges[node]) == 1:
            gene_net.edges[child] = gene_net.edges[node]
            del gene_net.edges[node]
            delIdx = gene_net.nodes.index(node)
            gene_net.nodes.pop(delIdx)
            gene_net.node_wts.pop(delIdx)
            num_networks -= 1

    # find leaves
    edges = {k: v for k, v in gene_net.edges.items()}
    internal_nodes = set([node for nodes in edges.values() for node in nodes])
    leaves = set(gene_net.edges) - internal_nodes

    # start DFS currently break ties randomly
    for leaf in leaves:
        curr = leaf
        while True:
            if curr in edges:
                # parents = edges[curr]
                parent = list(edges[curr])[0]
                del edges[curr]
                num_networks -= 1
                curr = parent
            else:
                break
    return num_networks


def rescue_UMI(eq_obj, min_set_txps, predictions):
    # first make new eqclass based on min_set_txps
    new_eqclasses = defaultdict(lambda: defaultdict(int))
    for i in range(eq_obj.num_eqclasses):
        try:
            labels = eq_obj.labels[i]
        except:
            print("ERROR: wrong Number of eqclasses provided")
            exit(1)
        umis = eq_obj.umis[i]

        new_label = []
        genes_set = set([])
        # add txp to new label only if in min set
        for label in labels:
            txp_name = eq_obj.txps_list[label]
            genes_set.add(eq_obj.txp_to_gene_dict[txp_name])
            if txp_name in min_set_txps:
                new_label.append(label)

        # Assign new label to new eqclass only if all txps were from one gene
        if len(genes_set) == 1:
            sorted_labels = tuple(sorted(new_label))
            for umi, count in umis.items():
                new_eqclasses[sorted_labels][umi] += count
        else:
            print("ERROR: {} Eqclass has txps from different gene", i)
            exit(1)
    print("New Alevin Eqclasses: ", new_eqclasses)

    # maintain a network of umis at each txp level for 1-length eqclass
    # Note No de-duplication yet
    txp_network = {}
    for labels, umis in new_eqclasses.items():
        if len(labels) == 1:
            txp_id = list(labels)[0]
            if txp_id in txp_network:
                print("ERROR: repetition of txp {} in new eqclasses".format(txp_id))
                exit(1)
            txp_network[txp_id] = Dedup.get_txp_umi_network(umis, txp_id)

    # making gene network
    gene_network = {}
    for gene in eq_obj.genes_list:
        txps = eq_obj.gene_to_txp_dict[gene]
        if gene in gene_network:
            print("ERROR: Duplicate gene network in Alevin_MST")
        net = Dedup.get_gene_umi_network(txps, txp_network)
        if net is not None:
            gene_network[gene] = net

    # Processing >1-length eqclasses
    for labels, umis in new_eqclasses.items():
        # Assuming all multi-g class are thrown away
        if len(labels) > 1:
            gene_name = eq_obj.txp_to_gene_dict[eq_obj.txps_list[labels[0]]]
            if gene_name in gene_network:
                Dedup.update_gene_network(umis, gene_network[gene_name], labels)
            else:
                # shouldn't be possible
                print("ERROR: SHOULDN'T REACH HERE (at least theoratically)")
                exit(1);

    for gene in eq_obj.genes_list:
        predictions[gene] = get_max_span_networks(gene_network[gene])


def check_set_coverage(eq_obj, txp_list):
    for eqclass in eq_obj.labels:
        covered = False
        for txp in eqclass:
            txp_name = eq_obj.txps_list[txp]
            if txp_name in txp_list:
                covered = True
                break
        if not covered:
            return False
    return True


def get_set_cover(eq_obj):
    # right_set is the set of txps
    right_set = eq_obj.txps_list

    # exhaustively search for all combination to get min_set cover
    grouping = 0
    min_txps = set([])
    while True:
        grouping = grouping + 1
        is_covered = False
        for txps in combinations(right_set, grouping):
            curr_min_txp_set = set(txps)
            is_covered = check_set_coverage(eq_obj, curr_min_txp_set)
            if is_covered:
                min_txps |= set(txps)
                break
        if is_covered:
            break
    print("min_set of txps for Alevin: ", min_txps)

    return min_txps


def get_prediction(eq_obj):
    predictions = {}
    for gene in eq_obj.genes_list:
        predictions[gene] = 0

    min_set_txps = get_set_cover(eq_obj)
    rescue_UMI(eq_obj, min_set_txps, predictions)

    return predictions
