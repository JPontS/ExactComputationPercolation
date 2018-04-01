import numpy as np

class CC:
    def __init__(self, nodes, roots, probability, sizes):
        self._nodes = np.sort(nodes)
        self._roots = np.array(sorted(set(roots)))
        self._CC_dict = {n: r for n,r in zip(nodes,roots)}
        self._probability = probability
        self._sizes = {r: s for r,s in zip(self._roots, sizes)}
        
    def clean_up(self):
        for node in self._nodes:
            self._CC_dict[node] = self.find_root(node)
        self._nodes = np.array(sorted(set(self._CC_dict.keys())))
        self._roots = np.array(sorted(set(self._CC_dict.values())))
        new_sizes = dict()
        for s_node in self._sizes.keys():
            new_sizes[self._CC_dict[s_node]] = new_sizes.get(self._CC_dict[s_node], np.zeros(101)) + self._sizes[s_node]
        self._sizes = new_sizes
        
    def find_root(self, node):
        return self._CC_dict[self._CC_dict[node]]
        
    def merge(self, CC2):
        self._sizes =  {root: self._sizes.get(root, np.zeros(101))*CC2._probability 
             + CC2._sizes.get(root, np.zeros(101))*self._probability for root in sorted(set(self._roots) | set(CC2._roots))}
        for node in CC2._nodes:
            if node in self._CC_dict:
                rs, rg = sorted([self.find_root(CC2._CC_dict[node]), self.find_root(node)])
                self._CC_dict[rg] = rs
            else:
                self._CC_dict[node] = CC2._CC_dict[node]
        self.clean_up()