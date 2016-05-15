import emf

#emf.run(r"G:\Projects\216041_Line111EMF\Models\FIELDS_inputs.xlsx",
#        r"G:\Projects\216041_Line111EMF\Models\new_model_output")

b = emf.run('practice_xcs.xlsx', path = 'run-dest/')

#b = emf.load_template('practice_xcs.xlsx')

#for xc in b:
#    xc.compare_DAT('XC-comparisons/' + xc.name.upper() + '.DAT',
#                    round = 3, path = 'XC-comparisons/')

#plt.show()
