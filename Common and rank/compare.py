import pandas as pd

# 1st version of comparing!!
# Works fine
name_1 = 'testA'  # -----------------
input_file_name_1 = 'BRAF monomer-dimer screens\\'+name_1+'.xlsx'
xls_1 = pd.ExcelFile(input_file_name_1)

print(xls_1.sheet_names)
print(len(xls_1.sheet_names))

name_2 = 'testB'  # --------------
input_file_name_2 = 'RAF1_ATPbiotin_Gavin\\'+name_2+'.xlsx'
xls_2 = pd.ExcelFile(input_file_name_2)

print(xls_2.sheet_names)
print(len(xls_2.sheet_names))

output_filename = 'comparison_Results\\' + name_1 + '_VS_' + name_2 + '.xlsx'
writer = pd.ExcelWriter(output_filename, engine='xlsxwriter')
for sheet_name_A in xls_1.sheet_names:
    #print(sheet_name_A)
    pd_A = xls_1.parse(sheet_name_A)
    pd_A.columns = [c.replace(' ', '_') for c in pd_A.columns]
    for sheet_name_B in xls_2.sheet_names:
        #print(sheet_name_B)
        pd_B = xls_2.parse(sheet_name_B)
        pd_B.columns = [c.replace(' ', '_') for c in pd_B.columns]
        # --------------------
        # ΠΡΟΣΟΧΗ ---Gene_names or Gene!!!!!
        common = pd_B[pd_B.Gene_names.isin(pd_A.Gene)]
        outpout_sheet_name = sheet_name_A + '__VS__' + sheet_name_B
        print(outpout_sheet_name, ' common: ', common.shape[0], ' from ', sheet_name_A, ':', pd_A.shape[0], ' and ', sheet_name_B, ':', pd_B.shape[0] )
        common.Gene_names.to_excel(writer, sheet_name=outpout_sheet_name[:31])
writer.close()
    
