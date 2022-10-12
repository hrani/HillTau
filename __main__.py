"""__main__.py: 
Entry point for this package.
"""
    
""" setup.py : Script for FindSim """
author__      = "HarshaRani"
__copyright__   = "Copyright 2018 FindSim, NCBS"
__maintainer__  = "HarshaRani"
__email__       = "hrani@ncbs.res.in"




def run():
    from hillhau.CppCode import hillTau
    hillTau.main()

def htgraph():
    from hilltau import htgraph
    htgraph.main()
 
def ht2sbml():
    from hilltau import ht2sbml
    ht2sbml.main()
 
if __name__ == '__main__':
    run()
    
