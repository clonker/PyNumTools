import os
import sys
import nose

import pynumtools

if __name__ == '__main__':
    target_dir = os.path.dirname(pynumtools.__file__)
    exit_status = nose.run(argv=[sys.argv[0], target_dir , '-v'])
    print("nose returned exit_status: ", exit_status)
    sys.exit(0 if exit_status else 1)
