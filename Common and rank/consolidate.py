import pandas as pd
import warnings
warnings.filterwarnings('ignore')


folderA = ["All interactors",
           "Differential interactors identified and quantified by log fold change",
           "\\searched for phosphosites\\Interactome dynamics of RAF1-BRAF_Fragpipe_phospho_lfq_summary"
           ]
folderB = ["RAF kinase assay mass spec results3 12-8-19",
           "RAF1 kinase assay_protein groups"
           ]

# =============================================================================
# folderA = [#"test_All interactors",
#            "test_Differential interactors identified and quantified by log fold change",
#            #"\\searched for phosphosites\\test_Interactome dynamics of RAF1-BRAF_Fragpipe_phospho_lfq_summary"
#            ]
# folderB = ["test_RAF kinase assay mass spec results3 12-8-19",
#            "test_RAF1 kinase assay_protein groups"
#            ]
# =============================================================================



#make big dataframe
big = pd.DataFrame(columns = ['Gene', 'proteinGroups', 'Raf1 Kinase specific', 'EF1-Raf1 Specific',
       'proteinGroups RAF kinase assay', 'RAF1', 'EF1-RAF1',
       'RAF1+BRAFV600E+DMSO', 'RAF1+BRAFV600E+Sor', 'RAF1+BRAFV600E+Vem',
       'RAF1+BRAFV600E+Dim', 'RAF1+BRAFV600E+Sor+Dim',
       'RAF1+BRAFV600E+Vem+Dim', 'RAF1+BRAF+DMSO', 'RAF1+BRAF+Sor',
       'RAF1+BRAF+Vem', 'RAF1+BRAF+Dim', 'RAF1+BRAF+Sor+Dim',
       'RAF1+BRAF+Vem+Dim', 'Sheet1'])
#big["Gene"] = ''
#BIG common file
big_filename = 'comparison_Results\\BIG.xlsx'
writer_big = pd.ExcelWriter(big_filename, engine='xlsxwriter')

# =============================================================================
# for name_1 in folderA:
#     input_file_name_1 = 'BRAF monomer-dimer screens\\'+name_1+'.xlsx'
#     xls_1 = pd.ExcelFile(input_file_name_1)
#     for name_2 in folderB:
#         input_file_name_2 = 'RAF1_ATPbiotin_Gavin\\'+name_2+'.xlsx'
#         xls_2 = pd.ExcelFile(input_file_name_2)
#         
#         for sheet_name_A in xls_1.sheet_names:
#             #print(sheet_name_A)
#             big[sheet_name_A] = ''
#             
#         for sheet_name_B in xls_2.sheet_names:
#             #print(sheet_name_B)
#             big[sheet_name_B] = ''
# =============================================================================
            
#print("cols ",big.columns)           

for name_1 in folderA:
    #name_1 = 'testA'  # -----------------
    input_file_name_1 = 'BRAF monomer-dimer screens\\'+name_1+'.xlsx'
    xls_1 = pd.ExcelFile(input_file_name_1)
    
    print("-----n1: ",name_1)
    #print(xls_1.sheet_names)
    #print(len(xls_1.sheet_names))
    for name_2 in folderB:
        print("-----n2: ",name_2)
        #name_2 = 'testB'  # --------------
        input_file_name_2 = 'RAF1_ATPbiotin_Gavin\\'+name_2+'.xlsx'
        xls_2 = pd.ExcelFile(input_file_name_2)
    
        #print(xls_2.sheet_names)
        #print(len(xls_2.sheet_names))
         
        if(name_1 == "\\searched for phosphosites\\Interactome dynamics of RAF1-BRAF_Fragpipe_phospho_lfq_summary"):
            output_filename = "comparison_Results\\Interactome dynamics"+ '_VS_' + name_2 + '.xlsx'
        else:
            output_filename = 'comparison_Results\\' + name_1 + '_VS_' + name_2 + '.xlsx'
        writer = pd.ExcelWriter(output_filename, engine='xlsxwriter')
        for sheet_name_A in xls_1.sheet_names:
            
            #print(sheet_name_A)
            pd_A = xls_1.parse(sheet_name_A)
            pd_A.columns = [c.replace(' ', '_') for c in pd_A.columns]
            
            if(name_1=="All interactors"):
                no_of_nulls = pd_A.Gene_names.isna().sum()
                pd_A = pd_A[pd_A.Gene_names.notnull()]
            else:
                no_of_nulls = pd_A.Gene.isna().sum()
                pd_A = pd_A[pd_A.Gene.notnull()]
            
            
            for sheet_name_B in xls_2.sheet_names:
                #print(sheet_name_B)
                pd_B = xls_2.parse(sheet_name_B)
                pd_B.columns = [c.replace(' ', '_') for c in pd_B.columns]
                #null delete
                no_of_nulls = pd_B.Gene_names.isna().sum()
                pd_B = pd_B[pd_B.Gene_names.notnull()]
                if(name_1=="All interactors"):
                    common = pd_B[pd_B.Gene_names.isin(pd_A.Gene_names)]
                else:
                    common = pd_B[pd_B.Gene_names.isin(pd_A.Gene)]
                
                outpout_sheet_name = sheet_name_A + '__VS__' + sheet_name_B
                #print(outpout_sheet_name, ' common: ', common.shape[0], ' from ', sheet_name_A, ':', pd_A.shape[0], ' and ', sheet_name_B, ':', pd_B.shape[0] )
                common.Gene_names.to_excel(writer, sheet_name=outpout_sheet_name[:31])
                #print("common:!!! ",common.Gene_names)
                for gene in common.Gene_names:
                    #print("geve ",gene)
                    xs = [None] * len(big.columns)
                    s = pd.Series(xs,index=big.columns)
                    if (big["Gene"] == gene).any():
                        l = big.index[big["Gene"]==gene].tolist()
                        #row = s[l]
                        #print("row= ",l)
                        big.iloc[l, big.columns.get_loc(sheet_name_B)] = "X"
                        big.iloc[l, big.columns.get_loc(sheet_name_A)] = "X"
                    else:
                        #print("NO ",gene)
                        s["Gene"] = gene
                        s[sheet_name_B] = "X"
                        s[sheet_name_A] = "X"
                        #df2 = pd.DataFrame(data=None, columns=big.columns)
                        #df2["Gene"] = gene
                        big = big.append(s,ignore_index=True)
                    
                
        writer.close()
big.to_excel(writer_big)
writer_big.close()
