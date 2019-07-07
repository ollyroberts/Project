import subprocess
print('Start Main onehelixres_refactored')
subprocess.call(['python', 'onehelixres_refactored.py', '../existing_code/1ct5.sec', '../existing_code/1ct5.1hr'])
print('Start Main twohelixres')
subprocess.call(['python', 'twohelixres_refactored.py', '../existing_code/1ct5.sec', '../existing_code/1h3l.2hr'])
print('Finish Main')