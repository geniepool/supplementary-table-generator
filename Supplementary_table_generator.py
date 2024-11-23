#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''

Requirements:

reference_genome file in GTF format - or any other human GTF, e.g. "hg38.refGene.gtf". Download links:
    hg38 - https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/
    hg19 - https://hgdownload.soe.ucsc.ed nu/goldenPath/hg19/bigZips/genes/
    T2T - https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/

Python libraries pandas, requests, pyarrow and fastparquet. install using:
pip3 install pandas requests pyarrow fastparquet

'''


import pandas as pd
import requests
import json

studied_gene = 'FABP5' # change to the identifier of your gene as listed in the GTF file. Notice - some GTF files do not use gene symbols (e.g. BRCA1), but other identifiers such as ENSEMBL IDs (e.g. ENSG00000012048)
referenceGenome = 'hg38' # The reference genome you are working with, e.g. hg19, hg38 or T2T
minCoverage = 20 # Minimal read coverage
minQual = 200 # Minimal quality score, taken from the VCF column "QUAL".
minAlphaMissenseScore = 0.564 # Missense variants with lower AlphaMissense score will be ommited

print(f'Working on {studied_gene}...')
gtf_file_path = 'hg38.refGene.gtf.gz' # change to your GTF file path

attributes_df = pd.read_parquet('attributes.parquet')
exons = pd.read_csv(gtf_file_path, header = None, sep = '\t')
exons['gene_id'] = exons[8].apply(lambda x : x.split('gene_id "')[1].split('"')[0])
exons = exons[exons['gene_id'] == studied_gene]
if exons.empty:
    print('studied_gene does not appear in the GTF file.')
    exit()
exons = exons[exons['gene_id'] == studied_gene]
exons = exons[exons[2] == 'exon']
exons = exons[[0,3,4]]
s3map = pd.read_csv('S3.map', sep = '\t', header = None) # include the exact path for files in the data lake
s3map.columns = ['version', 'reference', 'chromosome', 'modulus', 'path']

def listVariants(chromosome, position, mutations): # Digests the API's output and makes it suitable for a tabular output file
    lines = []
    coordinate = chromosome + ':' + str(position)
    for mutation in mutations:
        ref = mutation['ref']
        alt = mutation['alt']
        
        dbSNP = ''
        if 'dbSNP' in mutation.keys():
            if mutation['dbSNP'].startswith('rs'):
                dbSNP = mutation['dbSNP']
        
        AlphaMissense = ''
        if 'alphamissense' in mutation.keys():
            AlphaMissense = mutation['alphamissense']
        reference_coordinate = ''
        if f'{referenceGenome}_coordinate' in mutation.keys():
            reference_coordinate = mutation[f'{referenceGenome}_coordinate']
        
        variant = ref + '>' + alt
        homs = mutation['hom']
        hets = mutation['het']
        try:
            impact = mutation['impact']
        except:
            impact = ''
        lines.append([coordinate, variant, homs, hets, impact, dbSNP, AlphaMissense, reference_coordinate])
    return lines


def getAPI(coordinates): # Preformes the API request
    commonSamples = None
    coordinates = coordinates.upper().replace(' ', '').replace(',','').replace('CHR','').replace('MT','').strip()
    chromosome = coordinates.split(':')[0]
    pos = coordinates.split(':')[1]
    start, end = int(pos.split('-')[0]), int(pos.split('-')[1])
    try:
        query = f'http://geniepool-env-1.eba-ih62my9c.us-east-1.elasticbeanstalk.com/rest/index/{referenceGenome}/{coordinates}?ad={minCoverage}&qual={minQual}'
        data = requests.get(query).text 
        data = json.loads(data)
        if len(data) == 0:
            result = [html.P('No results')]
            return [result, n_clicks, {'display':'none'}, None, inputUpdate,'']
        elif len(data) == 1:
            lines = listVariants(chromosome, coordinates.split(':')[1].replace(',','').replace(' ',''), data['entries'])
            variantNumber = len(lines)
        else:
            variantNumber = int(data['count'])
            df = pd.json_normalize(data['data'])
            if df.empty:
                return None
    except:
        pos_range = range(start, end + 1)
        modulus1, modulus2 = start//100000, end//100000
        results = []
        for m in range(modulus1, modulus2 + 1):
            expression = "reference == @referenceGenome and chromosome == @chromosome and modulus == @m"
            query = s3map.query(expression)['path'].tolist()[0]
            query['S3'] = query.apply(lambda x : queryS3(x, pos_range), axis = 1)
            results += [i for i in query['S3'].tolist() if i.empty == False]
        if len(results) == 0:
            return None
        df = pd.concat(results)
        variantNumber = sum([len(i) for i in df['entries'].tolist()])
    data = df.apply(lambda x : listVariants(chromosome, x['pos'], x['entries']), axis = 1).tolist()
    lines = []
    for mutation in data:
        for line in mutation:
            hets = line[2]
            homs = line[3]
            if len(hets + homs) != 0:
                line[2] = hets
                line[3] = homs
                lines += [line]
    if len(lines) == 0:
            return None
        
    df = pd.DataFrame(lines, columns = ['Coordinate', 'Variant', 'Homozygote Samples', 'Heterozygote Samples', 'Impact', 'dbSNP', 'AlphaMissense', f'{referenceGenome}_coordinate'])
    del df[f'{referenceGenome}_coordinate']
    df['Homozygote Samples'] = df['Homozygote Samples'].apply(lambda x : [i['id'] for i in x if len(x) > 0])
    df['Heterozygote Samples'] = df['Heterozygote Samples'].apply(lambda x : [i['id'] for i in x if len(x) > 0])
    df = df[df['AlphaMissense'].apply(lambda x : '.' in str(x) or '1' in str(x))]
    df = df[df['AlphaMissense'] >= minAlphaMissenseScore]
    return df

queries = exons.apply(lambda x : x[0].replace('chr', '') + ':' + str(x[3]) + '-' + str(x[4]), axis = 1).tolist() # prepares the queries for the API
results = [getAPI(i) for i in queries]
df = pd.concat(results, axis=0) # Gathers all the results to a single dataframe
df.drop_duplicates(subset = ['Coordinate', 'Variant'], keep = 'first', inplace = True)
df['Homozygotes #'] = df['Homozygote Samples'].str.len()
df['Heterozygotes #'] = df['Heterozygote Samples'].str.len()

def getAttributes(relevant_samples): # Gets samples associated with a genomic variant and returns all associated attributes regarding them from BioSample
    attributes = attributes_df[attributes_df['Run'].isin(relevant_samples)]['Attributes'].tolist()
    attributes = sorted(set([item for sublist in attributes for item in sublist]))
    return attributes

df['Heterozygote attributes'] = df['Heterozygote Samples'].apply(getAttributes)
df['Homozygote attributes'] = df['Homozygote Samples'].apply(getAttributes)

# Output generation
df.to_csv(f'{studied_gene}.GeniePool.tsv', sep = '\t', index = False)
df.to_parquet(f'{studied_gene}.GeniePool.parquet', index = False)
print('Done')
