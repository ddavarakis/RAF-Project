import pandas as pd
import warnings
warnings.filterwarnings('ignore')


folderA = ["All interactors",
           #"Differential interactors identified and quantified by log fold change",
           "\\searched for phosphosites\\Interactome dynamics of RAF1-BRAF_Fragpipe_phospho_lfq_summary"
           ]

#BIG common file
big_filename = 'comparison_Results\\Rank_NEW.xlsx'
writer_big = pd.ExcelWriter(big_filename, engine='xlsxwriter')
      
big = pd.read_excel('comparison_Results\\BIG.xlsx')

rank = big.copy()
# insert new columns
rank.insert(1,"All-Rank","")
rank.insert(1,"Int-Rank","")
rank.insert(1,"Int-sumLFQ","")
rank.insert(1,"All-sumLFQ","")


rank = rank.drop('Unnamed: 0',axis=1)
rank = rank[['All-sumLFQ', 'Int-sumLFQ', 'Int-Rank', 'All-Rank', 'Gene',
       'proteinGroups', 'RAF1+BRAFV600E+DMSO', 'RAF1+BRAFV600E+Sor', 'RAF1+BRAFV600E+Vem',
       'RAF1+BRAFV600E+Dim', 'RAF1+BRAFV600E+Sor+Dim',
       'RAF1+BRAFV600E+Vem+Dim', 'RAF1+BRAF+DMSO', 'RAF1+BRAF+Sor',
       'RAF1+BRAF+Vem', 'RAF1+BRAF+Dim', 'RAF1+BRAF+Sor+Dim',
       'RAF1+BRAF+Vem+Dim', 'Sheet1', 'Raf1 Kinase specific', 'EF1-Raf1 Specific',
       'proteinGroups RAF kinase assay', 'RAF1', 'EF1-RAF1'
        ]]

print("rank columns: ", rank.columns)

for name_1 in folderA:
    #name_1 = 'testA'  # -----------------
    input_file_name_1 = 'BRAF monomer-dimer screens\\'+name_1+'.xlsx'
    xls_1 = pd.ExcelFile(input_file_name_1)
    print("-----n1: ",name_1)
    
    for sheet_name_A in xls_1.sheet_names:
        #print(sheet_name_A)
        pd_A = xls_1.parse(sheet_name_A)
        pd_A.columns = [c.replace(' ', '_') for c in pd_A.columns]
        
        df1 = pd_A.filter(regex='LFQ')
        pd_A["sumLFQ"] = df1.sum(axis = 1)
        for gene in big.Gene:        
            
            if(name_1 == "\\searched for phosphosites\\Interactome dynamics of RAF1-BRAF_Fragpipe_phospho_lfq_summary"):
                lnew = pd_A.index[pd_A["Gene"]==gene].tolist()
                colos = pd_A.columns.get_loc('sumLFQ')
                #print("col ",colos)
                #print("lnew = ",lnew," GENE ",gene)
                
                lfq = pd_A.iloc[lnew, colos]
                #print("lfq = ",lfq)
                
                lrank = rank.index[rank["Gene"]==gene].tolist()
                # Int-sumLFQ =========================== POSITION IN RANK =====================
                crank = rank.columns.get_loc('Int-sumLFQ')
                rank.iloc[lrank,crank] = lfq.max()
            else:
                lnew = pd_A.index[pd_A["Gene_names"]==gene].tolist()
                colos = pd_A.columns.get_loc('sumLFQ')
                #print("col ",colos)
                #print("lnew = ",lnew," GENE ",gene)
                
                lfq = pd_A.iloc[lnew, colos]
                #print("lfq = ",lfq)
                
                lrank = rank.index[rank["Gene"]==gene].tolist()
                # 'All-sumLFQ'=========================== POSITION IN RANK ===============
                crank = rank.columns.get_loc('All-sumLFQ')
                rank.iloc[lrank,crank] = lfq.max()
            
#print("",rank.columns)
rank.sort_values('Int-sumLFQ',inplace = True,ascending = False)
rank['Int-Rank'] = rank['Int-sumLFQ'].rank(ascending = False)
rank.sort_values('All-sumLFQ',inplace = True,ascending = False)
rank['All-Rank'] = rank['All-sumLFQ'].rank(ascending = False)

rank.rename(columns = {'Sheet1':'Interactome'}, inplace = True)

rank.to_excel(writer_big)
writer_big.close()
