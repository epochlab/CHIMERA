#!/usr/bin/env python3

class NodeTree(object):
    def __init__(self, left=None, right=None):
        self.left = left
        self.right = right

    def children(self):
        return self.left, self.right

def build_tree(nodes):
    while len(nodes) > 1:
        (key1, c1) = nodes[-1]
        (key2, c2) = nodes[-2]
        nodes = nodes[:-2]
        node = NodeTree(key1, key2)
        nodes.append((node, c1 + c2))
        nodes = sorted(nodes, key=lambda x: x[1], reverse=True)
    return nodes[0][0]

def assign_code(node, bin=''):
    if type(node) is str:
        return {node: bin}
    (l, r) = node.children()
    d = dict()
    d.update(assign_code(l, bin+'0'))
    d.update(assign_code(r, bin+'1'))
    return d

def encode(str, dict):
    output = ''
    for ch in str: output += dict[ch]
    return output