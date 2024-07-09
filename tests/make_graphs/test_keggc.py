import os

import networkx as nx

from kegg_pathways_completeness.bin.make_graphs.make_graphs \
    import recursive_parsing

dir_ = os.path.dirname(os.path.abspath(__file__))


def load_kegg_pathways():
    from importlib.resources import files
    data = files('kegg_pathways_completeness.pathways_data')
    c = data.joinpath('all_pathways.txt').open()
    cb = data.joinpath("all_pathways_names.txt").open()
    r = dict()
    b, cbb = c.readline(), cb.readline()
    while b:
        b, cbb = b.split(':'), cbb.split(':')
        assert b[0] not in r, b[0]
        r[b[0]] = (b[1][:-1], cbb[1][:-1])
        b, cbb = c.readline(), cb.readline()
    return r


def get_pathway_graph(pathway):
    """ Given KEGG pathway description return its graph representation
    constructed by the MGnify KEGG pathway completeness project
    """
    grh = nx.MultiDiGraph()
    grh.add_node(0)
    grh.add_node(1)
    grh, dict_edges, unnecessary_nodes = recursive_parsing(
        G=grh,
        dict_edges={},
        unnecessary_nodes=[],
        expression=pathway,
        start_node=0, end_node=1,
        weight=1)
    return grh


def is_pathway_complete(grh, abc):
    """ Given pathway graph, grh, and abundance values, abc,
        is the KEGG pathway represented by the graph complete? """
    grh_ = grh.copy()
    for c, r, b in list(grh.edges.data('label')):  # source, target, label
        # Delete edges when no annotations found
        if b not in abc and (c, r) in grh_.edges:
            grh_.remove_edge(c, r)

    # check whether we have routes connecting node-0 to node-1
    r = nx.has_path(grh_, 0, 1)
    nedges = len(grh_.edges)
    del grh_
    return r, nedges


def test_is_pathway_complete():
    pathways = load_kegg_pathways()
    abc = {
        "K00958": 123,
        "K00394": 23,
        "K11180": 10
    }
    pathway = get_pathway_graph(pathways['M00596'][0])
    c, _ = is_pathway_complete(pathway, abc)
    assert c is False
    abc = {
        "K00958": 123,
        "K00394": 23,
        "K00395": 23,
        "K11180": 10,
        "K11181": 10, "K27196": 10, "K27187": 8, "K27188": 8,
        "K27189":  8, "K27190":  8, "K27191": 5
    }
    c, _ = is_pathway_complete(pathway, abc)
    assert c is True
    pathway = get_pathway_graph(pathways['M00176'][0])
    c, _ = is_pathway_complete(pathway, abc)
    assert c is False
    # test with all required annotations missing
    c, nedges = is_pathway_complete(pathway, {})
    assert c is False
    assert nedges == 0
