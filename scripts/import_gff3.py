# -*- coding: utf-8 -*-

import pandas as pd
from collections import defaultdict

"""
Function that loads a gff file into a pandas dataframe
"""
def loadfile(filename,verbose):
        
    try:
        if verbose==True:
            print 'Loading', filename
        # obtaning sample names and number from 3rd line in header
        with open(filename) as f:
            f.readline()
            f.readline()
            sample_names = f.readline().split()[2].split(',')
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
    
        return data
    except:
        print 'Error loading the file'


"""
Function that check the header then load a gff file and check the content
Returns the dataframe if the format is ok, false if not
"""
def load_gff3(filename):
    try:
        Error = False     
        # Checking the format file 
        # Header and 1st data row
        header_rows=0
        with open(filename) as file:
            rowfile=file.readline()
            header_rows=header_rows+1

            while True:
                if rowfile.startswith('#'):
                    rowfile=file.readline()
                    header_rows=header_rows+1
                else:
                    data_1=rowfile.split('\t')

                    if header_rows!=4:
                        Error = True
                        print 'Bad header' 
                    break
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
                print('line', i, 'bad type error')
    
            # start<end
            if dataframe.loc[i, 'start'] >= dataframe.loc[i, 'end']:
                Error = True
                print('line', i, 'start >=end error')
    
            # Strand + or -
            if dataframe.loc[i, 'strand'] not in ['+', '-']:
                Error = True
                print('line', i, 'bad strand error')
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
        #Checking expression data
        expression_cols=[col for col in dataframe.columns if 'Expression_' in col]         
        for col in expression_cols:
            for i in range(dataframe.shape[0]):
                if not dataframe.loc[i,col].isdigit():
                    print dataframe.loc[i,col].isdigit()
                    print 'Expression count error line',i
                    Error= True
            dataframe[col]=dataframe[col].astype(int) #setting the datatype of counts
        if 'Filter' in dataframe.columns:         
            for i in range(dataframe.shape[0]):
                if dataframe.loc[i, 'Filter']!='Pass':
                    print 'Warning non-pass filter in line',i                   
        if Error:
            print 'File format error'
            return False
        print '--------------------------------------'
        print dataframe.dtypes
        print '--------------------------------------'
        print 'Format ok'
        return dataframe
    except:
        print 'Error checking the file'
        return False




file='/home/rafa/mirtop_project/Scripts/blood.gff3'
#file='/home/rafa/mirtop_project/Scripts/gff_examples/2uid_missing.gff'
#file='/home/rafa/mirtop_project/Scripts/gff_examples/coldata_missing.gff'
#file='/home/rafa/mirtop_project/Scripts/gff_examples/missing_filter_type.gff'
#file='/home/rafa/mirtop_project/Scripts/gff_examples/3wrong_type.gff'

filedata = load_gff3(file)



# Adquiring the sample names in a dataframe containing a gff3 file
def samplenames(dataframe):
    return [col for col in dataframe.columns if 'Expression_' in col]

# Determining the list of variants presente in a dataframe containing a gff3 file
def listvariants(filedata): 
    dataframe=filedata[filedata.type=='isomiR']
    list_variants_present=[]
    for index, row in dataframe.iterrows():      
        actual_variant=dataframe.loc[index,'Variant'].split(',')
        for i in range(len(actual_variant)):
            if actual_variant[i] not in list_variants_present:
                list_variants_present.append(actual_variant[i])
    return list_variants_present

# Number smallRNA by sample
def number_smallRNA_by_sample(dataframe):
    sample_names=samplenames(dataframe)
    print 'Number smallRNA by sample'
    print '--------------------------------------'
    total_number_dict={}
    for col in sample_names:
        total_number_dict[col]=dataframe[dataframe[col]>0][col].count()
        print  total_number_dict[col], col
    print '--------------------------------------'
    return total_number_dict

# Number of miRNA or isomiR in each sample
def number_smallRNA_by_type_by_sample(dataframe):
    sample_names=samplenames(dataframe)
    result=pd.DataFrame(index=sample_names, columns=['miRNA','isomiR'])
    for col in sample_names:
        result.at[col,'miRNA']=dataframe[(dataframe.type=='ref_miRNA')&(dataframe[col]>0)][col].count()
        result.at[col,'isomiR']=dataframe[(dataframe.type=='isomiR')&(dataframe[col]>0)][col].count()
    print result
    return result

# Number of isomir by type of modification
def number_isomir_by_type_by_sample(filedata):
    dataframe=filedata[filedata.type=='isomiR']
    sample_names=samplenames(dataframe)
    list_variants_present=listvariants(dataframe)
    result=pd.DataFrame(index=sample_names, columns=list_variants_present)
    result=result.fillna(0)
    for col in sample_names:
        for index,row in dataframe.iterrows():
            actual_variant=[]
            actual_variant=dataframe.loc[index,'Variant'].split(',')
            for i in range(len(actual_variant)):
                var=actual_variant[i]
                if dataframe.loc[index,col]!=0:
                    result.at[col,var]+=1
    print result
    return result


for col in samplenames(filedata):
    print filedata[filedata[col]>0][col].describe()

# Total counts by sample
def total_counts_by_sample(dataframe):
    sample_names=samplenames(dataframe)
    print 'Total counts by sample'
    print '--------------------------------------'
    total_counts_dict={}
    for col in sample_names:
        total_counts_dict[col]=dataframe[col].sum()
        print  total_counts_dict[col], col
    print '--------------------------------------'
    return total_counts_dict

# Counts by type miRNA or isomiR
def counts_by_type_by_sample(dataframe):
    sample_names=samplenames(dataframe)
    result=pd.DataFrame(index=sample_names, columns=['miRNA','isomiR'])
    for col in sample_names:
        result.at[col,'miRNA']=dataframe.groupby('type')[col].sum()['ref_miRNA']
        result.at[col,'isomiR']=dataframe.groupby('type')[col].sum()['isomiR']
    print result
    return result

# Counts by type of isomiR
def counts_by_isomir_type(filedata):
    dataframe=filedata[filedata.type=='isomiR']
    sample_names=samplenames(dataframe)
    list_variants_present=listvariants(dataframe)
    result=pd.DataFrame(index=sample_names, columns=list_variants_present)
    result=result.fillna(0)
    for col in sample_names:
        for index,row in dataframe.iterrows():
            actual_variant=[]
            actual_variant=dataframe.loc[index,'Variant'].split(',')
            for i in range(len(actual_variant)):
                var=actual_variant[i]
                result.at[col,var]+=dataframe.loc[index,col]
    print result
    return result



number_smallRNA_by_sample(filedata)
number_smallRNA_by_type_by_sample(filedata)
number_isomir_by_type_by_sample(filedata)
total_counts_by_sample(filedata)
counts_by_type_by_sample(filedata)
counts_by_isomir_type(filedata)

# Plots
import matplotlib.pyplot as plt
import numpy as np
list_variants_present=listvariants(filedata)
list_variants_present.sort()

# Counts of all variants in all samples
counts_type_isomir=counts_by_isomir_type(filedata)
counts_type_isomir=counts_type_isomir[list_variants_present]
fig=counts_type_isomir.T.plot(kind='line',title='Counts all modifications')
fig.set_xticks(np.arange(0,23,1))
fig.set_xticklabels(counts_type_isomir.columns, rotation=90)

# counts of 3' modifications in all samples
list_3=['iso_3p:+4','iso_3p:+3','iso_3p:+2','iso_3p:+1','iso_3p:-1','iso_3p:-2','iso_3p:-3','iso_3p:-4']
fig=counts_type_isomir[list_3].T.plot(kind='line',title='Counts with 3 end modifications')
fig.set_xticks(np.arange(0,8,1))
fig.set_xticklabels(list_3, rotation=90)

# counts of 5' modifications in all samples
list_5=['iso_5p:+3','iso_5p:+2','iso_5p:+1','iso_5p:-1','iso_5p:-2','iso_5p:-3']
fig=counts_type_isomir[list_5].T.plot(kind='line',title='Counts with 5 end modifications')
fig.set_xticks(np.arange(0,6,1))
fig.set_xticklabels(list_5, rotation=90)

# Number of isomir with all variants in all samples
num_isomir_type=number_isomir_by_type_by_sample(filedata)
num_isomir_type=num_isomir_type[list_variants_present]
fig=num_isomir_type.T.plot(kind='line',title='Number of IsomiR by all modifications')
fig.set_xticks(np.arange(0,23,1))
fig.set_xticklabels(num_isomir_type.columns, rotation=90)

# Number of isomir with 3' modifications in all samples
list_3=['iso_3p:+4','iso_3p:+3','iso_3p:+2','iso_3p:+1','iso_3p:-1','iso_3p:-2','iso_3p:-3','iso_3p:-4']
fig=num_isomir_type[list_3].T.plot(kind='line',title='Number of IsomiR with 3 end modifications')
fig.set_xticks(np.arange(0,8,1))
fig.set_xticklabels(list_3, rotation=90)

# Number of isomir with 5' modifications in all samples
list_5=['iso_5p:+3','iso_5p:+2','iso_5p:+1','iso_5p:-1','iso_5p:-2','iso_5p:-3']
fig=num_isomir_type[list_5].T.plot(kind='line',title='Number of IsomiR with 5 end modifications')
fig.set_xticks(np.arange(0,6,1))
fig.set_xticklabels(list_5, rotation=90)