import emf

#emf.run(r"G:\Projects\216041_Line111EMF\Models\FIELDS_inputs.xlsx",
#        r"G:\Projects\216041_Line111EMF\Models\new_model_output")

b = emf.load_template(r"G:\Projects\216041_Line111EMF\Models\FIELDS_inputs.xlsx")

emf.plot_groups(b, path = r"G:\Projects\216041_Line111EMF\Models\new_model_output\comparisons")

b = emf.load_template('practice_xcs.xlsx')

b['xc1'].compare_DAT(r"P:\MBaum\Programming\Python\python_code\FLD\XC-comparisons\XC1.DAT",
                    save = True, path = 'XC-comparisons/', round = 3)

#plt.show()
