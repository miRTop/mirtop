# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np


"""
Function that loads a gff file into a pandas dataframe
"""
def loadfile(filename,verbose,zero=False):
        
    try:
        if verbose==True:
            print 'Loading', filename
        # obtaning sample names and number from 3rd line in header
        with open(filename) as f:
            rowfile=f.readline()
            while True:
                if rowfile.startswith('## COLDATA'):
                    sample_names=rowfile.split()[2].split(',')
                    break
                else:
                    rowfile=f.readline()
        sample_number = len(sample_names)
        if verbose==True:
            print '--------------------------------------'
            print sample_number,' samples in the file'
            print '--------------------------------------'
            for elem in sample_names:
                print elem
            print '--------------------------------------'
        
        # number of columns in gff file
        gff_cols = pd.read_table(filename, sep='\t', skiprows=3, header=None).columns
         
        # Adquiring non-attributes data
        body_data=pd.read_table(filename, sep='\t', skiprows=3, header=None, usecols=gff_cols[0:-1])
        body_data.columns = ['SeqID', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase']
        # Adquiring attributes data
        atr_data = pd.read_table(filename, sep='\t', skiprows=3, header=None, usecols=gff_cols[[-1]])
    
        # Splitting the attributes column
        list_atr = []
        # cheking attributes present in first row
        attr_names=atr_data.values[0,0].split()[0::2]      #attributes in the column
        num_attr = len(attr_names)                         #number of attributes
        expression_colindex=attr_names.index ('Expression')  #position of the expression column in the attr column 
        if verbose==True:
            print num_attr,' attributes in the file '
            print '--------------------------------------'
            for attr in attr_names:
                print attr
            print '--------------------------------------'
        # joining rows of attributes without the descriptor
        for row in range(atr_data.shape[0]):
            list_atr.append(atr_data.values[row, 0].split()[1::2])
    
        # appending observations
        atr_data = pd.DataFrame(list_atr, columns=attr_names)
        # Eliminating ;
        for i in attr_names:
            atr_data[i] = atr_data[i].replace({';': ''}, regex=True)
    
        # desglosing the expression column in a column for each sample
        list_expression=[]
        for row in range(atr_data.shape[0]):
            list_expression.append(atr_data.values[row,expression_colindex].split(','))
        sample_names=['Expression_' +x for x in sample_names]
        
        expression_data = pd.DataFrame(list_expression, columns=sample_names)      
        
        
        atr_data=atr_data.drop('Expression',axis=1) #Remove the expression column
        atr_data=atr_data.join(expression_data)
    
        # Joining the body and attributes dataframes
        data = body_data.join(atr_data)
        
        #Unlisting the variant column
        tempframe=data[data.type=='isomiR']
        list_variants_present=[]
        for row in tempframe.itertuples():      
            actual_variant=tempframe.loc[row.Index,'Variant'].split(',')
            for i in range(len(actual_variant)):
                if actual_variant[i] not in list_variants_present:
                    list_variants_present.append(actual_variant[i])
                    
                    
        
        for var in list_variants_present:
            for row in data.itertuples():
                if var in data.loc[row.Index,'Variant']:
                    #print var, data.loc[row.Index,'Variant']
                    data.at[row.Index,var]=1
                else:
                    data.at[row.Index,var]=np.nan    
        return data
    except:
        print 'Error loading the file'


"""
Function that check the header then load a gff file and check the content
Returns the dataframe if the format is ok, false if not
"""
def load_check_gff3(filename):
    try:
        Error = False
        coldata_found=False
        # Checking the format file 
        # Header and 1st data row
        with open(filename) as file:
            rowfile=file.readline()           
            while True:
                if rowfile.startswith('##'):
                    if rowfile.startswith("## COLDATA"):
                         coldata_found=True
                    rowfile=file.readline()
                else:
                    data_1=rowfile.split('\t')
                    break
        if coldata_found==False:
            print 'No COLDATA, bad header'
            return False
        
        #Number of columns without breaking down attributes column 
        if len(data_1) > 9:
            Error=True
            print(len(data_1))
            print('Too much columns')
        # Cheking the attributes column  
        attr_names=data_1[-1].split(';')
        list_attr=[]
        for atr in range(len(attr_names)-1):
            list_attr.append(attr_names[atr].split()[0])
        possible_attr = ['UID', 'Read', 'Name', 'Parent', 'Variant', 'Cigar', 'Hits', 'Alias', 'Genomic', 'Expression',
                           'Filter', 'Seed_fam']
        for attr in list_attr:    
            if attr not in possible_attr:
                Error=True
                print attr,'is not a possible attribute'
                break  
        if Error:
            print 'File format error'            
            return False      
        # If not format error, loading content
        try:
            dataframe=loadfile(filename,True)
        except:
            print 'Error loading file'
            return False
        print 'Checking content'
        for i in range(dataframe.shape[0]):
            # Labels in type column
            if dataframe.loc[i, 'type'] not in ['ref_miRNA', 'isomiR']:
                Error = True
                print'line', i, 'bad type error'
    
            # start<end
            if dataframe.loc[i, 'start'] >= dataframe.loc[i, 'end']:
                Error = True
                print 'line', i, 'start >=end error'
    
            # Strand + or -
            if dataframe.loc[i, 'strand'] not in ['+', '-']:
                Error = True
                print 'line', i, 'bad strand error'
            # Variant checking
            possible_variant=['iso_5p','iso_3p','iso_add','iso_snp_seed','iso_snp_central_offset','iso_snp_central',
                  'iso_central_supp','iso_snp_central_supp','iso_snp']
            
            variant_i=dataframe.loc[i,'Variant'].split(',')
            
            if len(variant_i)==1 and variant_i[0]!='NA':
                if variant_i[0].split(':')[0] not in possible_variant:
                    Error = True
                    print 'Variant error', variant_i[0].split(':')[0], 'line', i
            elif variant_i[0]!='NA':
                for var in range(len(variant_i)):
                    if variant_i[var].split(':')[0] not in possible_variant:
                        Error = True
                        print 'Variant error', variant_i[0].split(':')[0], 'line', i
        #dataframe=dataframe.drop('Variant',index=1) #Remove the Variant column
        #Checking expression data
        expression_cols=[col for col in dataframe.columns if 'Expression_' in col]         
        for col in expression_cols:
            for i in range(dataframe.shape[0]):
                if not dataframe.loc[i,col].isdigit():
                    print dataframe.loc[i,col].isdigit()
                    print 'Expression count error line',i
                    Error= True
            dataframe[col]=dataframe[col].astype(int) #setting the datatype of counts
            dataframe[col]=dataframe[col].replace(0,np.nan) #Setting 0 reads to NaN
        if 'Filter' in dataframe.columns:         
            for i in range(dataframe.shape[0]):
                if dataframe.loc[i, 'Filter']!='Pass':
                    print 'Warning non-pass filter in line',i                   
        if Error:
            print 'File format error'
            return False
        dataframe=dataframe.drop('Variant',axis=1) #Remove the Variant column
        print '--------------------------------------'
        print dataframe.dtypes
        print '--------------------------------------'
        print 'Format ok'
        return dataframe
    except:
        print 'Error checking the file'
        return False




#file='/home/rafa/mirtop_project/Scripts/blood.gff3'
#file='/home/rafa/mirtop_project/Scripts/blood2.gff3'
file='/home/rafa/mirtop_project/Aligments/blood/blood2/13_plasma-1-TruSeq_fastq-mirbase-ready.gff3'

#file='/home/rafa/mirtop_project/Scripts/gff_examples/2uid_missing.gff'
#file='/home/rafa/mirtop_project/Scripts/gff_examples/coldata_missing.gff'
#file='/home/rafa/mirtop_project/Scripts/gff_examples/missing_filter_type.gff'
#file='/home/rafa/mirtop_project/Scripts/gff_examples/3wrong_type.gff'

filedata = load_check_gff3(file)
