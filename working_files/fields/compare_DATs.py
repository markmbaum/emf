import os
import emf.fields as fld

sb = fld.load_template('practice_xcs.xlsx')

for xs in sb:
    dat_path = os.path.join('compare_DAT', xs.sheet.upper() + '.DAT')
    if(os.path.isfile(dat_path)):
        xs.compare_DAT(dat_path, path='compare_DAT')
