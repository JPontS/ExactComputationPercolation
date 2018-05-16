import numpy as np

def log_factorial(x):    
    return np.sum(np.log(np.arange(1,x+1,1)))

def compute_Binomial(N, n):
    return log_factorial(N) - log_factorial(N-n) - log_factorial(n)

def find_root_general(roots, node):
    while node != roots[node]:
        node=roots[node]
    return node

class CC:
    def __init__(self, nodes, roots, probability, sizes):
        self._nodes = np.sort(nodes)
        self._tuple_roots = roots
        self._roots = np.array(sorted(set(roots)))
        self._CC_dict = {n: r for n,r in zip(nodes,roots)}
        self._probability = probability
        self._sizes = {r: s for r,s in zip(self._roots, sizes)}
        
    def clean_up(self):
        for node in self._nodes:
            self._CC_dict[node] = self.find_root(node)
        self._tuple_roots = tuple(self._CC_dict[node] for node in self._nodes)
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
        self._probability = self._probability*CC2._probability
        self._nodes = np.array(sorted(set(self._nodes) | set(CC2._nodes)))
        self.clean_up()
        
    def remove_node_from_CC(self, node):
        rooted_nodes = sorted({ k for k,v in self._CC_dict.items() if v == node })
        if not rooted_nodes:
            self._sizes[self._CC_dict[node]] = self._sizes[self._CC_dict[node]] + self._probability
        if len(rooted_nodes) > 1:
            new_root = rooted_nodes[1]
            for item in rooted_nodes[1:]:
                self._CC_dict[item] = new_root
            self._sizes[new_root] = self._sizes[node] + self._probability
        self._sizes.pop(node, None)
        self._CC_dict.pop(node, None)   
        self._nodes = sorted(self._CC_dict.keys())
        self._tuple_roots = tuple(self._CC_dict[node] for node in self._nodes)
        self._roots = np.array(sorted(set(self._tuple_roots)))
        
class CCD:
    def __init__(self, nodes, ls_roots):
        self._nodes = np.sort(nodes)
        self._roots = {roots: np.array(sorted(set(roots))) for roots in ls_roots}
        self._probabilities = {roots: np.abs(i-np.arange(0.0, 1.01, 0.01)) for i,roots in enumerate(ls_roots)}
        self._sizes = {roots: np.zeros((self._roots[roots].shape[0], 101)) for roots in ls_roots} 
        
    def multiply(self, CCD2):
        probabilities, sizes, roots = dict(), dict(), dict()
        k1, k2 = list(self._probabilities.keys()), list(CCD2._probabilities.keys())
        for CC1_roots in k1:
            for CC2_roots in k2:
                CC1 = CC(self._nodes, CC1_roots, self._probabilities[CC1_roots], self._sizes[CC1_roots])
                CC2 = CC(CCD2._nodes, CC2_roots, CCD2._probabilities[CC2_roots], CCD2._sizes[CC2_roots])
                CC1.merge(CC2)
                probabilities[CC1._tuple_roots] = probabilities.get(CC1._tuple_roots, np.zeros(101)) + CC1._probability
                sizes[CC1._tuple_roots] = sizes.get(CC1._tuple_roots, np.zeros((CC1._roots.shape[0], 101))) \
                                            + np.vstack([value for (key, value) in sorted(CC1._sizes.items())])
                roots[CC1._tuple_roots] = CC1._roots
        self._probabilities, self._sizes, self._roots = probabilities, sizes, roots
        self._nodes = np.array(sorted(set(self._nodes) | set(CCD2._nodes)))
        
    def remove_node(self, node):
        probabilities, sizes, roots = dict(), dict(), dict()
        for CC_roots in self._probabilities.keys():
            CC_new = CC(self._nodes, CC_roots, self._probabilities[CC_roots], self._sizes[CC_roots]) 
            CC_new.remove_node_from_CC(node)
            probabilities[CC_new._tuple_roots] = probabilities.get(CC_new._tuple_roots, np.zeros(101)) + CC_new._probability 
            sizes[CC_new._tuple_roots] = sizes.get(CC_new._tuple_roots, np.zeros((CC_new._roots.shape[0], 101))) \
                                                + np.vstack([value for (key, value) in sorted(CC_new._sizes.items())])
            roots[CC_new._tuple_roots] = CC_new._roots
        self._probabilities, self._sizes, self._roots = probabilities, sizes, roots        
        self._nodes = np.array(sorted(set(self._nodes) - {node}))  

    def delete_CCD(self):
        self._nodes = np.zeros(1,int)
        self._roots = np.zeros(1, int)
        self._probabilities, self._sizes = {}, {}

def Newman_Ziff_algorithm(edges, N, perc_type, node_i=None, node_j=None):
    
    N_nodes, N_edges = np.max(edges), len(edges)
    Q, Q_mean  = np.zeros((N_edges, N)), np.zeros(N_edges + 1)
        
    for k in range(N):
        np.random.shuffle(edges)
        roots = np.arange(1, N_nodes+1, 1, int)
        Size, roots = dict(zip(roots, np.ones((N_nodes), int))), dict(zip(roots, roots))

        for n,bond in enumerate(edges):
            root1 , root2 = find_root_general(roots, bond[0]), find_root_general(roots, bond[1])
            if root1 != root2 :
                rg, rs = (root1, root2) if Size[root1] > Size[root2] else (root2, root1)
                roots[rs] = rg
                Size[rg], Size[rs] = Size[root1] + Size[root2], 0
            
            if perc_type == 'pair_wise_i_j' and find_root_general(roots, node_i) == find_root_general(roots, node_j):
                Q[n,k]= 1

            elif perc_type == 'size_node_i':
                Q[n,k]= Size[find_root_general(roots,node_i)]

            elif perc_type == 'total_size':
                Q[n,k]= sum(x*x for x in Size.values())/len(Size.values())
    
    Q_mean[0] = 1.0 if perc_type != 'pair_wise_i_j' else 0.0
    Q_mean[1:] = np.mean(Q, axis=1)

    B_n = np.zeros(N_edges+1)
    for n in range(1,N_edges+1):
        B_n[n]=compute_Binomial(N_edges,n)

    p_range = np.arange(0.0,1.01,0.01)
    Q_p = np.zeros(p_range.shape)
    Q_p[0] = 1.0 if perc_type != 'pair_wise_i_j' else 0.0
    
    for i,p in enumerate(p_range[1:-1],1):
        for j,element in enumerate(Q_mean):
            Q_p[i] = Q_p[i] + element*np.exp(B_n[j] + j*np.log(p) + (N_edges - j)*np.log(1.0 - p))

    Q_p[-1] = N_nodes if perc_type != 'pair_wise_i_j' else 1.0
    return Q_p