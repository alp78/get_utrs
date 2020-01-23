import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import mygene
from gtfparse import read_gtf
import pybedtools
import logging
from os import getcwd, path, remove
from sys import exit
from termcolor import colored

def make_utr_tables(gtf, logger):
    logger.info(colored('Extracting UTRs from GTF...', 'blue'))
    gtf_df = read_gtf(gtf)
    # 3'UTR
    utr3_df = gtf_df[gtf_df['feature']=='three_prime_utr']
    utr3_unique = utr3_df.drop_duplicates(subset=['gene_id', 'source'], keep='first')
    utr3_unique.reset_index(drop=True, inplace=True)
    utr3_unique = utr3_unique.loc[:, ['source', 'gene_id', 'feature', 'seqname', 'strand', 'start', 'end']]
    utr3_unique.to_csv('utr3.txt', sep='\t', index = False)
    # 5'UTR
    utr5_df = gtf_df[gtf_df['feature']=='five_prime_utr']
    utr5_unique = utr5_df.drop_duplicates(subset=['gene_id', 'source'], keep='first')
    utr5_unique.reset_index(drop=True, inplace=True)
    utr5_unique = utr5_unique.loc[:, ['source', 'gene_id', 'feature', 'seqname', 'strand', 'start', 'end']]
    utr5_unique.to_csv('utr5.txt', sep='\t', index = False)

def convert_id(_id, logger):
    logger.info(colored('Retrieving Ensembl ID...', 'blue'))
    mg = mygene.MyGeneInfo()
    entrez_id = mg.query(_id, species='human')['hits'][0]['entrezgene']
    gene_symbol = mg.query(entrez_id, species='human')['hits'][0]['symbol']
    ensembl_id = mg.getgene(entrez_id, 'ensembl')['ensembl']['gene']
    return ensembl_id, gene_symbol

def get_3utr(ensembl_id, gene_symbol, logger):
    logger.info(colored('Getting 3-UTR...', 'blue'))
    utr3_unique = pd.read_csv('utr3.txt', sep='\t', skiprows=(0), header=(0), low_memory=False)
    choices = utr3_unique[utr3_unique['gene_id']==ensembl_id]['source'].to_list()
    if len(choices) > 0:
        temp_df = utr3_unique[utr3_unique['gene_id']==ensembl_id]
        if 'ensembl' in choices:
            location = temp_df[temp_df['source']=='ensembl'].values.tolist()[0][-6:]
        elif 'ensembl_havana' in choices:
            location = temp_df[temp_df['source']=='ensembl_havana'].values.tolist()[0][-6:]
        else:
            location = temp_df[temp_df['source']=='havana'].values.tolist()[0][-6:]
        location[2] = 'chr' + location[2]
        location[0] = gene_symbol
    return location

def get_5utr(ensembl_id, gene_symbol, logger):
    logger.info(colored('Getting 5-UTR...', 'blue'))
    utr5_unique = pd.read_csv('utr5.txt', sep='\t', skiprows=(0), header=(0), low_memory=False)
    choices = utr5_unique[utr5_unique['gene_id']==ensembl_id]['source'].to_list()
    if len(choices) > 0:
        temp_df = utr5_unique[utr5_unique['gene_id']==ensembl_id]
        if 'ensembl' in choices:
            location = temp_df[temp_df['source']=='ensembl'].values.tolist()[0][-6:]
        elif 'ensembl_havana' in choices:
            location = temp_df[temp_df['source']=='ensembl_havana'].values.tolist()[0][-6:]
        else:
            location = temp_df[temp_df['source']=='havana'].values.tolist()[0][-6:]
        location[2] = 'chr' + location[2]
        location[0] = gene_symbol
    return location

def make_utrs_fasta(ref_fasta_full_path, location_list, logger):
    logger.info(colored('Cooking FASTA file...', 'blue'))
    recs_list = []
    for location in location_list:
        print(location)
        strand = location[3]
        loc = pybedtools.BedTool(f"""{location[2]} {location[4]} {location[5]}""", from_string=True)
        fasta = pybedtools.example_filename(ref_fasta_full_path)
        seq = loc.sequence(fi=fasta)
        out_fasta = seq.save_seqs('temp_UTR.fasta')
        for rec in SeqIO.parse('temp_UTR.fasta', 'fasta'):
            seq = rec.seq.upper()
            if strand == '-':
                rev_seq = seq.reverse_complement()
                rec.seq = rev_seq
            rec.id = f'{location[0]}|{location[1]}|{location[2]}|{location[3]}|{location[4]}|{location[5]}'
            rec.description = ''
            recs_list.append(rec)
    recs_gen = (rec for rec in recs_list)
    SeqIO.write(recs_gen, 'UTR.fasta', "fasta")
    remove('temp_UTR.fasta')

# PROGRAM START
def get_utrs(gtf, ref_fasta_full_path, _id):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    location_list = []
    ensembl_id = ''
    gene_symbol = ''
    utr3 = []
    utr5 = []

    if not (path.isfile('utr3.txt') and path.isfile('utr5.txt')):
        try:
            make_utr_tables(gtf, logger)
        except:
            logger.info(colored('Error in parsing GTF file', 'red'))
            del logger
            exit(0)

    try:
        ensembl_id, gene_symbol = convert_id(_id, logger)
    except:
        logger.info(colored('Error in retrieving Ensembl ID', 'red'))
        del logger
        exit(0)

    if ensembl_id:
        try:
            utr3 = get_3utr(ensembl_id, gene_symbol, logger)
        except:
            logger.info(colored('Error in retrieving 3-UTR', 'red'))
            del logger
            exit(0)
        try:   
            utr5 = get_5utr(ensembl_id, gene_symbol, logger)
        except:
            logger.info(colored('Error in retrieving 5-UTR', 'red'))
            del logger
            exit(0)  
    else:
        logger.info(colored('Error in retrieving Ensembl ID', 'red'))
        del logger
        exit(0)

    location_list.append(utr5)
    location_list.append(utr3)

    if (len(location_list) == 2 and len(utr3) == 6 and len(utr5) == 6):
        try:
            make_utrs_fasta(ref_fasta_full_path, location_list, logger)
        except:
            logger.info(colored('Error in creating fasta file', 'red'))
            del logger
            exit(0)          
    else:
        logger.info(colored('Error in location list', 'red'))
        del logger
        exit(0)

    logger.info(colored('All done!', 'green'))
    del logger
