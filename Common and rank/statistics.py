import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

#files = ["BRAF monomer-dimer screens\\All interactors"]
files = ["RAF1_ATPbiotin_Gavin\\RAF1 kinase assay_protein groups"]

for name_1 in files:
    #name_1 = 'All interactors'  # -----------------
    input_file_name_1 = name_1+'.xlsx'
    xls_1 = pd.ExcelFile(input_file_name_1)

    
    print("------- Filename: ", name_1, "----------")


    for sheet_name_A in xls_1.sheet_names:
        #print(sheet_name_A)
        print("---Sheet: ", sheet_name_A)

        lfq1 = pd.DataFrame()
        
        pd_A = xls_1.parse(sheet_name_A)
        print("Number of rows: ", pd_A.shape[0])

        print("Number of columns: ", pd_A.shape[1])

        pd_A.columns = [c.replace(' ', '_') for c in pd_A.columns]
        # check Gene column
        no_of_nulls = pd_A.Gene_names.isna().sum()
        print("Number of null values: ", no_of_nulls)

        pd_A = pd_A[pd_A.Gene_names.notnull()]
            
        df1 = pd_A.filter(regex='LFQ')
        lfq1["Gene"] = pd_A.Gene_names
        lfq1["sumLFQ"] = df1.sum(axis = 1)
        
        
        a_duplicates = lfq1[lfq1['Gene'].duplicated(keep=False) == True]
        print("Total Number of duplicate values: ", a_duplicates.shape[0])

        dupl = lfq1.pivot_table(columns=['Gene'], aggfunc='size')
        dupl.to_frame()
        
        #uniquedup = dupl[dupl>1][0]
        ld = dupl[dupl>1].index
        
        print("Unique duplicates are: ", ld)

        a_duplicates = a_duplicates.sort_values('Gene')
        #print("SORTED duplicates are: ", a_duplicates)
    







    
