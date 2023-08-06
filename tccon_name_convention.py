## Imported libraries ##
import os
import re
import shutil
## Definition of data path ##

opusfiles = os.listdir('/home/STG05_TCCON/work/test_i2s/altzomoni/2022/')

## Change names to TCCON convention standars ##

for opusfile in opusfiles:
	
	if re.search(r'\w*.\d+',opusfile):
		name = re.search(r'(\w+)SN.(\d+)',opusfile)
		
		date = name.group(1)
		number = str(name.group(2).zfill(3))
		print(date)
		new_name  = 'al20{}s0e00a.{}'.format(date,number)
		print(new_name)
		shutil.copy(os.path.join('/home/STG05_TCCON/work/test_i2s/altzomoni/2022/', opusfile),new_name)