import networkx as nx
import random as rand
import connectome_utils as utl


def algorithm_wrapper(args):
    """

    :param args: algorithm (fn), graph (networkx), filepath (str, for output)
    :return:
    """
    algorithm, graph, filepath = args[:3]
    if len(args) > 3:
        kwargs = args[3]
    else:
        kwargs = dict()

    result = algorithm(graph, **kwargs)

    utl.json_serialise(result, filepath)
    return 0


def di_maslov(graph, swap_prop=10, return_successful=False):
    """
    Based on matlab code here http://www.cmth.bnl.gov/~maslov/dir_generate_srand.m

    :param graph: graph to be randomised (not in place)
    :param swap_prop: How many swaps should be attempted per edge
    :param return_successful: Return a tuple including what proportion of attempted edge swaps were actually completed
    :return: randomised, binarised graph
    """

    count = 0
    edgeset = set(graph.edges_iter())
    for _ in range(len(edgeset)*swap_prop):
        ((src1, tgt1), (src2, tgt2)) = rand.sample(edgeset, 2)
        if not {(src1, tgt2), (src2, tgt1)} & (edgeset - {(src1, tgt1), (src2, tgt2)}):
            edgeset.remove((src1, tgt1))
            edgeset.remove((src2, tgt2))
            edgeset.add((src1, tgt2))
            edgeset.add((src2, tgt1))
            count += 1

    out_graph = nx.DiGraph()
    out_graph.add_nodes_from(graph.nodes_iter())
    out_graph.add_edges_from(edgeset)

    if return_successful:
        return out_graph, count / (len(edgeset) * swap_prop)
    else:
        return out_graph


def undi_maslov(graph, swap_prop=10, return_successful=False):
    """
    Based on matlab code here http://www.cmth.bnl.gov/~maslov/dir_generate_srand.m

    :param graph: graph to be randomised (not in place)
    :param swap_prop: How many swaps should be attempted per edge
    :param return_successful: Return a tuple including what proportion of attempted edge swaps were actually completed
    :return: randomised, binarised graph
    """

    count = 0
    edgeset = set(nx.Graph(graph).edges_iter())
    for _ in range(len(edgeset)*swap_prop):
        ((src1, tgt1), (src2, tgt2)) = rand.sample(edgeset, 2)
        if not {(src1, tgt2), (src2, tgt1), (tgt2, src1), (tgt1, src2)} &\
                (edgeset - {(src1, tgt1), (tgt1, src1), (src2, tgt2), (tgt2, src2)}):
            edgeset.remove((src1, tgt1))
            edgeset.remove((src2, tgt2))
            edgeset.add((src1, tgt2))
            edgeset.add((src2, tgt1))
            count += 1

    out_graph = nx.Graph()
    out_graph.add_nodes_from(graph.nodes_iter())
    out_graph.add_edges_from(edgeset)

    if return_successful:
        return out_graph, count / (len(edgeset) * swap_prop)
    else:
        return out_graph