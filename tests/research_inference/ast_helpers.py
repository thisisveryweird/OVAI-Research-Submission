"""Offline extraction: no source imports, imports/decorators/top-level calls executed."""
import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

def source(name):
    return (ROOT / 'scripts' / name).read_text(encoding='utf-8-sig')

def extract(name, functions, globals=None):
    tree = ast.parse(source(name))
    nodes = [ast.ImportFrom(module='__future__', names=[ast.alias(name='annotations')], level=0)]
    constants = {}
    def safe_constant(value):
        if isinstance(value, ast.Name) and value.id in constants:
            return constants[value.id]
        if isinstance(value, ast.BinOp) and isinstance(value.op, ast.Add):
            left, right = safe_constant(value.left), safe_constant(value.right)
            if isinstance(left, tuple) and isinstance(right, tuple):
                return left + right
            raise ValueError('Only tuple concatenation allowed')
        return ast.literal_eval(value)
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            value = node.value
            try:
                constant = safe_constant(value)
            except (ValueError, TypeError):
                continue
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if not all(isinstance(target, ast.Name) for target in targets):
                continue
            for target in targets:
                constants[target.id] = constant
            nodes.append(node)
        elif isinstance(node, ast.FunctionDef) and node.name in functions:
            node.decorator_list = []
            nodes.append(node)
    found = {n.name for n in nodes if isinstance(n, ast.FunctionDef)}
    if set(functions) - found:
        raise AssertionError('Missing contract helpers: ' + str(set(functions) - found))
    namespace = dict(globals or {})
    exec(compile(ast.fix_missing_locations(ast.Module(body=nodes, type_ignores=[])), name, 'exec'), namespace)
    return namespace
