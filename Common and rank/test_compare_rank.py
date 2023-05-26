import pandas as pd
import warnings
warnings.filterwarnings('ignore')



folder = ["All interactors_VS_RAF kinase assay mass spec results3 12-8-19",
          "All interactors_VS_RAF1 kinase assay_protein groups",
          "Differential interactors identified and quantified by log fold change_VS_RAF kinase assay mass spec results3 12-8-19",
          "Differential interactors identified and quantified by log fold change_VS_RAF1 kinase assay_protein groups", 
          "Interactome dynamics_VS_RAF kinase assay mass spec results3 12-8-19",
          "Interactome dynamics_VS_RAF1 kinase assay_protein groups"
           ]


big = pd.read_excel('comparison_Results\\BIG.xlsx')

# Check                                        
for name in folder:
        input_file_name = 'comparison_Results\\'+name+'.xlsx'
        xls = pd.ExcelFile(input_file_name)
        print("-----n1: ",name)
        files = name.split('_VS_')
        
        for sheet_name in xls.sheet_names:
            print(sheet_name)
            columns = sheet_name.split('__VS__')
            #print(columns)
            pd_comp = pd.read_excel(xls, sheet_name)
  
            if (columns[1].startswith("RA")):
                column_2 = "RAF1"
                #print("RAF1")
          
            if (columns[1].startswith("Ra")):
                column_2 = "Raf1 Kinase specific"
                #print("Raf1 Kinase specific")
                

                
            if columns[1].startswith("EF") and files[1].startswith("RAF kinase assay mass spec"):
                column_2 = "EF1-Raf1 Specific" 
                #print("EF1-Raf1 Specific")
                
            if columns[1].startswith("EF") and files[1].startswith("RAF1 kinase assay_protein"):
                column_2 = "EF1-RAF1"
                #print("EF1-RAF1")

                
            if (columns[1].startswith("pro")):
                #print("PPPPPPPPPPPPPP")
                column_2 = "proteinGroups RAF kinase assay"    

            #print("columns: ", columns)
            
            for gene in pd_comp["Gene_names"]:
                if (big["Gene"] == gene).any():
                    l = big.index[big["Gene"]==gene].tolist()
                    if big.iloc[l][columns[0]].any() == 'X' :
                        big.iloc[l, big.columns.get_loc(columns[0])] = "OK"
                        #print("XXXXXXXXXX")
                    #print("Column:----", column_2, "----")
                    if big.iloc[l][column_2].any() == 'X':
                        big.iloc[l, big.columns.get_loc(column_2)] = "OK"
                else:
                    print("GENE - NOT FOUND????????????????????????: ", gene)                

big_filename = 'comparison_Results\\BIG-TESTED.xlsx'
writer_big = pd.ExcelWriter(big_filename, engine='xlsxwriter')
big.to_excel(writer_big)
writer_big.close()                              
            
            