"""Read-only dependency readiness probe; run with -B. No installation."""
import importlib
import json
import sys
import traceback


def main():
    print(json.dumps({'executable': sys.executable, 'version': sys.version}))
    failures = []
    for name in ('streamlit', 'joblib', 'numpy', 'pandas', 'openpyxl', 'sklearn'):
        try:
            module = importlib.import_module(name)
            print(json.dumps({'module': name, 'status': 'ok', 'version': getattr(module, '__version__', None)}))
        except Exception:
            failures.append(name)
            print(json.dumps({'module': name, 'status': 'failed'}))
            traceback.print_exc()
    return bool(failures)


if __name__ == '__main__':
    sys.exit(main())
