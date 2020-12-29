import subprocess
from subprocess import Popen, PIPE
import sys


def main():
    with open(sys.argv[3], 'r') as supposed:
        supposed_out = supposed.read().strip().split()

    args = [sys.executable, sys.argv[1], sys.argv[2]]
    with Popen(args, stdout=PIPE) as proc:
        try:
            out, err = proc.communicate(None, timeout=1)
            out = str(out, encoding='utf-8')
        except subprocess.TimeoutExpired:
            print('Too slow')
            proc.kill()
            return
        out = out.strip().split()

    for i, (a, b) in enumerate(zip(out, supposed_out)):
        if a != b:
            print(f'Result #{i} is wrong: {a} != {b}')
            return

    print('OK')


if __name__ == '__main__':
    main()
