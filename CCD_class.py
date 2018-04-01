import numpy as np

class CC:
    def __init__(self, nodes, roots):
        self._nodes = np.sort(nodes)
        self._roots = roots 
        self._CC_dict = {n: r for n,r in zip(nodes,roots)} 
        
    def clean_up(self):
        for node in self._nodes:
            self._CC_dict[node] = find_root(self._CC_dict,node)
    
    def find_roots(self, node):
        return self._CC_dict[self._CC_dict[node]]
    
    def merge(self, nodes2, CC2):
        for node in nodes2:
            if node in self._CC_dict:
                rs, rg = sorted([find_root(self._CC_dict,CC2[node]),find_root(self._CC_dict,node)])
                self._CC_dict[rg] = rs
            else:
                self._CC_dict[node] = CC2[node]
        self.clean_up()
        self._nodes = np.sort(list(set(nodes) | set(nodes2)))
    