from collections import defaultdict

#change this function based on the use - case
# currently assuming u3 is 1-edit of u1 and u4
def one_distance(u1, u2):
    first = int(u1.replace("u", ""))
    second = int(u2.replace("u", ""))
    if abs(first-second) <= 1:
        return True
    else:
        return False

class Network:
    def __init__(self, umis, tid):
        self.nodes = list([uid+"_"+str(tid) for uid in umis.keys()])
        self.edges = defaultdict(set)
        self.node_wts = list(umis.values())
        self.populate_edges(umis, tid)

    def populate_edges(self, umis, tid):
        visit_list = sorted(umis, key=umis.get)

        while(len(visit_list) > 0):
            removed_umis = [ visit_list.pop() ]

            while len(removed_umis) != 0:
                query_umi = removed_umis.pop(0)
                query_freq = umis[query_umi]

                for umi in visit_list[::-1]:
                    umi_freq = umis[umi]
                    if one_distance(umi, query_umi)  and ((query_freq/2.0)+1)>umi_freq:
                        removed_umis.append(umi)

                for umi in removed_umis:
                    self.edges[umi+"_"+str(tid)].add( query_umi+"_"+str(tid) )
                    visit_list.remove(umi)
        print ("Edges for {} in Alevin_MST".format(tid))
        print (self.edges)

    def combine_nets(self, other_net):
        self.nodes += other_net.nodes
        self.node_wts += other_net.node_wts
        for k,v in other_net.edges.items():
            if k in self.edges:
                print ("ALEVIN MST: network can't be combined")
            self.edges[k] = v

    def one_distance(self, first_txp, second, txps):
        toks = first_txp.split("_")
        if len(toks) == 1:
            return one_distance(first_txp, second)
        else:
            return one_distance(toks[0], second)

def get_txp_umi_network(umis, txp):
    return Network(umis, txp)

def get_gene_umi_network(txps, txp_net):
    net = None
    for txp in txps:
        if txp in txp_net:
            other_net = txp_net[txp]
            if net is None:
                net = other_net
            else:
                net.combine_nets(other_net)
    return net

def update_gene_network(umis, gene_net, txps):
    nodes = gene_net.nodes
    # traverse over each umi in multi-t eqclass
    for umi, count in umis.items():
        # check if the umi is in 1-edit from any node in the graph
        for nid, node in enumerate(nodes):
            # if within 1-edit then add relevant edge
            if gene_net.one_distance(node, umi, set(txps)):
                node_freq = gene_net.node_wts[nid]
                if (node_freq/2.0)+1 > count:
                    gene_net.edges[umi].add(node)
                elif (count/2.0)+1 > node_freq:
                    gene_net.edges[node].add(umi)
    gene_net.nodes += umis.keys()
    gene_net.node_wts += umis.values()

def collapse(visit_list, umis, pcr_correct):
    removed_umis = [ visit_list.pop() ]
    num_collapsed = 0

    while len(removed_umis) != 0:
        query_umi = removed_umis.pop(0)
        query_freq = umis[query_umi]

        for umi in visit_list[::-1]:
            umi_freq = umis[umi]
            if one_distance(umi, query_umi):
                if not pcr_correct or \
                        ( pcr_correct and ((query_freq/2.0)+1)>umi_freq ):
                    removed_umis.append(umi)
                    num_collapsed += 1

        for umi in removed_umis:
            visit_list.remove(umi)

    return num_collapsed

def dedup(umis, pcr_correct):
    '''
    deduplicate based on the umi and it's frequency
    :param umis: dictionary of the UMI and it's count
    :return: counts after deduplicating
    '''
    sorted_umis_seqs = sorted(umis, key=umis.get)
    num_molecules = len(sorted_umis_seqs)

    while(len(sorted_umis_seqs) > 0):
        num_collapsed = collapse(sorted_umis_seqs, umis, pcr_correct)
        num_molecules = num_molecules - num_collapsed
        if num_molecules < 1:
            print("ERROR in UMI collapsing")

    return num_molecules
