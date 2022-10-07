from Bio import Phylo
from io import StringIO

def tree_to_str(t):
    t_str = StringIO()
    Phylo.draw_ascii(t, t_str)
    t_str.seek(0)
    t_str.read()

def trees_approx_equal(t1_path: str, t2_path: str):
    t1 = Phylo.read(t1_path, 'newick')
    t1.ladderize()
    t2 = Phylo.read(t2_path, 'newick')
    t2.ladderize()
    return tree_to_str(t1) == tree_to_str(t2)
