from dataclasses import make_dataclass
import typing
import dataclasses
import json
import polyfempy

def getter(contain, key):
    if key not in contain:
        return []
    return contain[key]

def parse_tree(rules, mark_required = False):
    tree = dict()
    for r in rules:
        pt = r['pointer']
        pt = pt.replace('/*','')

        path = pt.split('/')
        if pt == '/':
            path = ['']

        subtree = tree
        for k in path:
            if k not in subtree:
                subtree[k] = dict()
            subtree = subtree[k]
        req, opt = (getter(r,'required'), getter(r, 'optional'))
        for k in req:
            if k not in subtree:
                subtree[k] = {}
            if mark_required:
                subtree[k]['R'] = '0'
        for k in opt:
            if k not in subtree:
                subtree[k] = {}
    return tree


def snake_to_camel(name):
    if name == '':
        return 'PolyFEM'
    return ''.join([n.capitalize() for n in name.split('_')])

def print_script():
    # deprecated to the make_dataclasses interface
    key = ''
    node = tree[key]
    indent = 4
    repre = 'from dataclasses import dataclass\nimport typing\n'

    def rec(key, val, h):
        global repre
        if key == 'lambda':
            return
        if key !='':
            repre += (' '*h*indent + f'{key} : typing.Any = None\n')

        if len(val) == 0:
            return
        classname = snake_to_camel(key)
        repre += (' '*h*indent + f'@dataclass\n')
        repre += (' '*h*indent + f'class {classname}:\n')
        for k,v in (val.items()):
            rec(k,v, h+1)

    rec(key, node, 0)

    repre += """
    if __name__ == '__main__':
        print(PolyFEM(geometry='mesh.obj'))
    """

    with open('test.py', 'w') as fp:
        fp.write(repre)



def recursive_data_class_builder(node, name):
    fields = []
    namespace = dict()
    for k, v in node.items():
        if k == 'lambda':
            continue
        fields.append((k,typing.Any,None))
        if isinstance(v, dict) and len(v) != 0:
            name0 = snake_to_camel(k)
            namespace[name0] = recursive_data_class_builder(v, name0)
    return make_dataclass(name, fields=fields, namespace = namespace)


def prune_none_entry(node):
    if not isinstance(node, dict):
        if isinstance(node, list):
            for v in node:
                prune_none_entry(v)
    else:
        for k, v in list(node.items()):
            prune_none_entry(v)
            if v is None or (hasattr(v, '__len__') and len(v) == 0):
                node.pop(k)

def generate_input_dict(config):
    output = dataclasses.asdict(config)
    prune_none_entry(output)
    return output


with open('data/default_rules.json')as fp:
    rules = json.load(fp)

tree = parse_tree(rules)
pf = recursive_data_class_builder(tree[''], '')

geometry = pf.Geometry(mesh='import.obj')
matnu = pf.Materials(id = 0, nu=0.1)
matE = pf.Materials(id = 1,E=10)

config = pf(geometry = geometry, materials = [matnu, matE])

def solve(**kwargs):
    d = generate_input_dict(config)
    polyfempy.solve(**d)

if __name__ == '__main__':
  