import os
from subprocess import Popen

import grip
import nbconvert
from odpydoc import doc

#generate documentation html files
doc('emf', script='ga.js')

#remove the old index file
if(os.path.isfile('index.html')):
    os.remove('index.html')
    print('index.html removed')

#make emf.html the index file
os.rename('emf.html', 'index.html')
print('emf.html -> index.html')

#render the readmes for fields and subcalc
grip.api.export(path='README-fields.md', out_filename='README-fields.html')
grip.api.export(path='README-subcalc.md', out_filename='README-subcalc.html')

#render the main readme
grip.api.export(path='../README.md', out_filename='../README.html')

#convert notebooks to html by calling batch file
print('Exporting notebooks to html')
p = Popen('autonotebookexport.bat')
stdout, stderr = p.communicate()
