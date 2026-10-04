"""Pure reusable static evidence; never imports or executes reviewed code."""
import ast
import hashlib

def source_evidence(content: bytes):
    text = content.decode('utf-8-sig')
    tree = ast.parse(text)
    return {'sha256': hashlib.sha256(content).hexdigest(), 'bytes': len(content),
            'lines': len(text.splitlines()),
            'functions': [{'name': n.name, 'start': n.lineno, 'end': n.end_lineno}
                          for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]}
