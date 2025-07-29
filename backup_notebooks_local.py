import os
import os.path as op
import shutil
from datetime import datetime
from subprocess import call
## sync from leffe
print ('Copying from leffe...')
src = f'buzzanga@leffe:/home/buzzanga/BB_LIB/VLM/bzFRInGE/notebooks/'
dst = f'./notebooks_leffe'
cmd = f'rsync -avzP {src} {dst}'
print (cmd)
call(cmd.split())


print ('\nBacking up to home for druva...')
today = datetime.today().strftime('%Y%m%d')
srcd  = dst
dstd  = op.join(op.expanduser('~'), 'backup', 'FRInGE_notebooks', today)

## use rsync cuz of symlinks
cmd = f'rsync -avzP {srcd} {dstd}'
print (cmd)
call(cmd.split())

