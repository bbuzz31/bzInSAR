import os
import os.path as op
import shutil
from datetime import datetime

today = datetime.today().strftime('%Y%m%d')
srcd  = op.join(os.getcwd(), 'notebooks')
dstd  = op.join(op.expanduser('~'), 'notebooks_backup', today)

shutil.copytree(srcd, dstd, dirs_exist_ok=True)


