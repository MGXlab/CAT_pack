#!/usr/bin/env python3

import argparse
import datetime
import gzip
import multiprocessing
import os
import subprocess
import sys
import tarfile
import urllib.request

import about
import check
import shared
import tax


def parse_arguments():
    date = str(datetime.datetime.now().date())
    
    parser = argparse.ArgumentParser(prog='CAT prepare',
                                     description='Download and construct '
                                                 'CAT/BAT database files.',
                                     usage='CAT prepare (--fresh | '
                                           '--existing) [options] [-h / '
                                           '--help]',
                                     add_help=False)
    
    required_choice = parser.add_argument_group('Required choice')

    group = required_choice.add_mutually_exclusive_group(required=True)

    group.add_argument('--fresh',
                       dest='fresh',
                       action='store_true',
                       help='Start with a fresh database.')
    group.add_argument('--existing',
                       dest='fresh',
                       action='store_false',
                       help='Start with an existing database. CAT prepare '
                            'will search the supplied database and taxonomy '
                            'folders and only construct files that do not '
                            'exist yet.')
    
    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument('-d',
                          '--database_folder',
                          dest='database_folder',
                          metavar='',
                          required=False,
                          type=str,
                          default='{0}_CAT_database'.format(date),
                          help='Name of folder to which database files will '
                               'be written (default: {date}_CAT_database).')
    optional.add_argument('-t',
                          '--taxonomy_folder',
                          dest='taxonomy_folder',
                          metavar='',
                          required=False,
                          type=str,
                          default='{0}_taxonomy'.format(date),
                          help='Name of folder to which taxonomy files will '
                               'be downloaded (default: {date}_taxonomy).')
    optional.add_argument('--path_to_diamond',
                          dest='path_to_diamond',
                          metavar='',
                          required=False,
                          type=str,
                          default='diamond',
                          help='Path to DIAMOND binaries. Please supply if CAT'
                               ' prepare can not find DIAMOND.')
    optional.add_argument('-q',
                          '--quiet',
                          dest='quiet',
                          required=False,
                          action='store_true',
                          help='Suppress verbosity.')
    optional.add_argument('--no_log',
                          dest='no_log',
                          required=False,
                          action='store_true',
                          help='Suppress log file.')
    optional.add_argument('-h',
                          '--help',
                          action='help',
                          help='Show this help message and exit.')
    
    specific = parser.add_argument_group('DIAMOND specific optional arguments')

    specific.add_argument('-n',
                          '--nproc',
                          dest='nproc',
                          metavar='',
                          required=False,
                          type=int,
                          default=multiprocessing.cpu_count(),
                          help='Number of cores to deploy by DIAMOND makedb '
                               '(default: maximum).')
    
    (args, extra_args) = parser.parse_known_args()

    extra_args = [arg for (i, arg) in enumerate(extra_args) if
                  (i, arg) != (0, 'prepare')]
    if len(extra_args) > 0:
        sys.exit('error: too much arguments supplied:\n{0}'
                 ''.format('\n'.join(extra_args)))

    return (args, date)


def download_taxonomy_files(taxonomy_folder, date, log_file, quiet):
    url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
    tmp_taxonomy_file = '{0}/{1}.taxdump.tar.gz'.format(taxonomy_folder, date)

    message = ('Downloading and extracting taxonomy files from {0} to {1}.'
               ''.format(url, taxonomy_folder))
    shared.give_user_feedback(message, log_file, quiet)
    
    try:
        urllib.request.urlretrieve(url, tmp_taxonomy_file)
    except:
        message = 'ERROR: donwload of taxonomy files failed.'
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)
    
    try:
        with tarfile.open(tmp_taxonomy_file) as tar:
            tar.extractall(taxonomy_folder)
    except:
        message = ('ERROR: something went wrong while extracting the taxonomy '
                   'files.')
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)
        
    message = 'Download complete!'
    shared.give_user_feedback(message, log_file, quiet)
    
    
def download_prot_accession2taxid_file(prot_accession2taxid_file,
                                       date,
                                       log_file,
                                       quiet):
    url = ('ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/'
           'prot.accession2taxid.gz')

    message = ('Downloading mapping file from {0} to {1}.'
               ''.format(url, prot_accession2taxid_file))
    shared.give_user_feedback(message, log_file, quiet)
    
    try:
        urllib.request.urlretrieve(url, prot_accession2taxid_file)
    except:
        message = 'ERROR: download of prot.accession2taxid.gz failed.'
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    message = 'Download complete!'
    shared.give_user_feedback(message, log_file, quiet)
    
    return prot_accession2taxid_file


def download_nr(nr_file, log_file, quiet):
    url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz'
    
    message = ('Downloading nr database from {0} to {1}.'
               ''.format(url, nr_file))
    shared.give_user_feedback(message, log_file, quiet)
    
    try:
        urllib.request.urlretrieve(url, nr_file)
    except:
        message = 'ERROR: download of nr database failed.'
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    message = 'Download complete!'
    shared.give_user_feedback(message, log_file, quiet)


def make_diamond_database(path_to_diamond,
                          nr_file,
                          diamond_database_prefix,
                          nproc,
                          log_file,
                          quiet):
    message = ('Constructing DIAMOND database {0}.dmnd from {1}. '
               'Please be patient...'.format(diamond_database_prefix, nr_file))
    shared.give_user_feedback(message, log_file, quiet)

    command = [path_to_diamond, 'makedb',
               '--in', nr_file,
               '-d', diamond_database_prefix,
               '-p', str(nproc),
               '--quiet']
    try:
        subprocess.check_call(command)
    except:
        message = 'ERROR: DIAMOND database could not be created.'
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)
        
    message = 'DIAMOND database constructed!'
    shared.give_user_feedback(message, log_file, quiet)
    
    
def import_prot_accession2taxid(prot_accession2taxid_file, log_file, quiet):
    message = ('Loading {0} into memory. Please be patient...'
               ''.format(prot_accession2taxid_file))
    shared.give_user_feedback(message, log_file, quiet)
    
    prot_accession2taxid = {}

    with gzip.open(prot_accession2taxid_file, 'rb') as f1:
        for line in f1:
            line = line.decode('utf-8')
            line = line.split('\t')

            prot_accession2taxid[line[1]] = line[2]

    return prot_accession2taxid


def make_fastaid2LCAtaxid_file(taxonomy_folder,
                               fastaid2LCAtaxid_file,
                               nr_file,
                               prot_accession2taxid_file,
                               log_file,
                               quiet):
    prot_accession2taxid = import_prot_accession2taxid(prot_accession2taxid_file,
                                                       log_file,
                                                       quiet)
    nodes_dmp = '{0}/nodes.dmp'.format(taxonomy_folder)
    (taxid2parent, taxid2rank) = tax.import_nodes(nodes_dmp, log_file, quiet)

    message = ('Finding LCA of all protein accession numbers in fasta headers '
               'of {0}. Please be patient...'.format(nr_file))
    shared.give_user_feedback(message, log_file, quiet)
    
    corrected = 0
    total = 0
    with gzip.open(nr_file, 'rb') as f1, open(fastaid2LCAtaxid_file, 'w') as outf1:
        for line in f1:
            line = line.decode('utf-8')
            if not line.startswith('>'):
                continue

            line = line.lstrip('>').split('\x01')

            accession_numbers = [i.split(' ')[0] for i in line]
            fastaid = accession_numbers[0]
            
            list_of_lineages = []
            for accession_number in accession_numbers:
                try:
                    taxid = prot_accession2taxid[accession_number]
                    lineage = tax.find_lineage(taxid, taxid2parent)
                    list_of_lineages.append(lineage)
                except:
                    # This accounts for missing accession numbers in
                    # prot.accession2taxid and missing nodes in nodes.dmp.
                    continue

            total += 1
            
            if len(list_of_lineages) == 0:
                # This accounts for entries that only contain accession numbers
                # that are missing in prot.accession2taxid or whose taxid is
                # missing in nodes.dmp. Note that these entries are thus not
                # present in the output file.
                continue

            LCAtaxid = tax.find_LCA(list_of_lineages)

            outf1.write('{0}\t{1}\n'.format(fastaid, LCAtaxid))

            try:
                if LCAtaxid != prot_accession2taxid[fastaid]:
                    corrected += 1
            except:
                # If the fastaid cannot be found in prot.accession2taxid, but
                # a taxid is given to the fastaid based on secondary accession
                # numbers, it is counted as a correction as well.
                corrected += 1

    message = ('Done! File {0} is created. '
               '{1} of {2} headers ({3:.1f}%) corrected. Please wait '
               'patiently for Python to collect carbage.'
               ''.format(fastaid2LCAtaxid_file,
                         corrected,
                         total,
                         corrected / total * 100))
    shared.give_user_feedback(message, log_file, quiet)
    
    
def find_offspring(taxonomy_folder, fastaid2LCAtaxid_file, log_file, quiet):
    nodes_dmp = '{0}/nodes.dmp'.format(taxonomy_folder)
    (taxid2parent, taxid2rank) = tax.import_nodes(nodes_dmp, log_file, quiet)

    message = 'Searching nr database for taxids with multiple offspring.'
    shared.give_user_feedback(message, log_file, quiet)

    taxid2offspring = {}

    with open(fastaid2LCAtaxid_file, 'r') as f1:
        for line in f1:
            line = line.rstrip().split('\t')

            taxid = line[1]
            lineage = tax.find_lineage(taxid, taxid2parent)

            for (i, taxid) in enumerate(lineage):
                # The first taxid in the lineage does not have a daughter node.
                if i == 0:
                    continue

                if taxid not in taxid2offspring:
                    taxid2offspring[taxid] = set()

                offspring = lineage[i - 1]

                taxid2offspring[taxid].add(offspring)
                
    return taxid2offspring


def write_taxids_with_multiple_offspring_file(taxids_with_multiple_offspring_file,
                                              taxid2offspring,
                                              log_file,
                                              quiet):
    message = 'Writing {0}.'.format(taxids_with_multiple_offspring_file)
    shared.give_user_feedback(message, log_file, quiet)

    with open(taxids_with_multiple_offspring_file, 'w') as outf1:
        for taxid in taxid2offspring:
            if len(taxid2offspring[taxid]) >= 2:
                outf1.write('{0}\n'.format(taxid))
                        
                        
def prepare(step_list,
            taxonomy_folder,
            database_folder,
            date,
            prot_accession2taxid_file,
            nr_file,
            path_to_diamond,
            diamond_database_prefix,
            nproc,
            fastaid2LCAtaxid_file,
            taxids_with_multiple_offspring_file,
            log_file,
            quiet):
    if 'download_taxonomy_files' in step_list:
        download_taxonomy_files(taxonomy_folder, date, log_file, quiet)

    if 'download_prot_accession2taxid_file' in step_list:
        download_prot_accession2taxid_file(prot_accession2taxid_file,
                                           date,
                                           log_file,quiet)
        
    if 'download_nr' in step_list:
        download_nr(nr_file, log_file, quiet)

    if 'make_diamond_database' in step_list:
        make_diamond_database(path_to_diamond,
                              nr_file,
                              diamond_database_prefix,
                              nproc,
                              log_file,
                              quiet)

    if 'make_fastaid2LCAtaxid_file' in step_list:
        make_fastaid2LCAtaxid_file(taxonomy_folder,
                                   fastaid2LCAtaxid_file,
                                   nr_file,
                                   prot_accession2taxid_file,
                                   log_file,
                                   quiet)

    if 'make_taxids_with_multiple_offspring_file' in step_list:
        taxid2offspring = find_offspring(taxonomy_folder,
                                     fastaid2LCAtaxid_file,
                                     log_file,
                                     quiet)
        write_taxids_with_multiple_offspring_file(taxids_with_multiple_offspring_file,
                                                  taxid2offspring,
                                                  log_file,
                                                  quiet)

    message = ('\n-----------------\n\n'
               '[{0}] CAT prepare is done!'.format(datetime.datetime.now()))
    shared.give_user_feedback(message, log_file, quiet, show_time=False)

    if nr_file is not None:
        message = 'You may remove {0} now.'.format(nr_file)
        shared.give_user_feedback(message, log_file, quiet, show_time=False)

    message = ('\nSupply the following arguments to CAT or BAT if you want to '
               'use the constructed database:\n'
               '-d / --database_folder {0}\n'
               '-t / --taxonomy_folder {1}'
               ''.format(database_folder,
                         taxonomy_folder))
    shared.give_user_feedback(message, log_file, quiet, show_time=False)
    
    
def run_fresh(args, date):
    (database_folder,
     taxonomy_folder,
     path_to_diamond,
     quiet,
     no_log,
     nproc) = check.convert_arguments(args)
    
    if no_log:
        log_file = None
    else:
        log_file = '{0}.CAT_prepare.fresh.log'.format(date)
        with open(log_file, 'w') as outf1:
            pass
        
    message = '# CAT v{0}.'.format(about.__version__)
    shared.give_user_feedback(message, log_file, quiet, show_time=False)
    
    message = ('\n'
               'CAT prepare is running, constructing a fresh database.\n'
               'Rawr!\n\n'
               'WARNING: preparing the database files may take a couple of '
               'hours.\n\n'
               'Supplied command: {0}\n\n'
               'Taxonomy folder: {1}/\n'
               'Database folder: {2}/\n'
               'Log file: {3}\n\n'
               '-----------------\n'.format(' '.join(sys.argv),
                                            taxonomy_folder,
                                            database_folder,
                                            log_file))
    shared.give_user_feedback(message, log_file, quiet, show_time=False)
    
    # Check diamond path.
    error = check.check_diamond_binaries(path_to_diamond, log_file, quiet)
    if error:
        sys.exit(1)
    
    # Check taxonomy folder.
    taxonomy_folder_inspect = check.inspect_taxonomy_folder(taxonomy_folder)
    if taxonomy_folder_inspect != [None]:
        if len([file for file in taxonomy_folder_inspect if
                file is not None]) > 0:
            message = ('ERROR: taxonomy folder {0} exists already and '
                       'contains taxonomy files. Please supply a novel or '
                       'empty folder if you want to start fresh, or run '
                       'CAT prepare --existing.'
                       ''.format(taxonomy_folder))
            shared.give_user_feedback(message, log_file, quiet, error=True)

            sys.exit(1)

        message = ('Taxonomy folder exists already. Taxonomy files will be '
                   'downloaded to it.')
        shared.give_user_feedback(message, log_file, quiet)
            
    database_folder_inspect = check.inspect_database_folder(database_folder)

    # Check database folder.
    if database_folder_inspect != [None]:
        if len([file_ for file_ in database_folder_inspect if
                file_ is not None]) > 0:
            message = ('ERROR: database folder {0} exists already and '
                       'contains database files. Please supply a novel or '
                       'empty folder if you want to start fresh.'
                       ''.format(database_folder))
            shared.give_user_feedback(message, log_file, quiet, error=True)

            sys.exit(1)

        message = ('Database folder exists already. Database file will be '
                   'downloaded to it / constructed in it.')
        shared.give_user_feedback(message, log_file, quiet)
        
    # Check memory.
    min_mem = 100
    (total_memory, error) = check.check_memory(min_mem)

    if error:
        message = ('ERROR: at least {0}GB of memory is needed for a fresh '
                   'database construction. {1}GB is found on your system. You '
                   'can either try to find a machine with more memory, or '
                   'download preconstructed database files from '
                   'tbb.bio.uu.nl/bastiaan/CAT_prepare/.'
                   ''.format(min_mem, total_memory))
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)
        
    if taxonomy_folder_inspect == [None]:
        os.mkdir(taxonomy_folder)

        message = '{0} is created.'.format(taxonomy_folder)
        shared.give_user_feedback(message, log_file, quiet)
    if database_folder_inspect == [None]:
        os.mkdir(database_folder)

        message = '{0} is created.'.format(database_folder)
        shared.give_user_feedback(message, log_file, quiet)
        
    prot_accession2taxid_file = ('{0}/{1}.prot.accession2taxid.gz'
                                 ''.format(taxonomy_folder, date))
    nr_file = '{0}/{1}.nr.gz'.format(database_folder, date)
    diamond_database_prefix = '{0}/{1}.nr'.format(database_folder, date)
    fastaid2LCAtaxid_file = ('{0}/{1}.nr.fastaid2LCAtaxid'
                             ''.format(database_folder, date))
    taxids_with_multiple_offspring_file = ('{0}/{1}.nr.taxids_with_multiple_'
                                           'offspring'
                                           ''.format(database_folder, date))
    
    step_list = ['download_taxonomy_files',
                 'download_prot_accession2taxid_file',
                 'download_nr',
                 'make_diamond_database',
                 'make_fastaid2LCAtaxid_file',
                 'make_taxids_with_multiple_offspring_file']
    
    prepare(step_list,
            taxonomy_folder,
            database_folder,
            date,
            prot_accession2taxid_file,
            nr_file,
            path_to_diamond,
            diamond_database_prefix,
            nproc,
            fastaid2LCAtaxid_file,
            taxids_with_multiple_offspring_file,
            log_file,
            quiet)
    
    
def run_existing(args, date):
    step_list = []

    (database_folder,
     taxonomy_folder,
     path_to_diamond,
     quiet,
     no_log,
     nproc) = check.convert_arguments(args)
    
    if no_log:
        log_file = None
    else:
        log_file = '{0}.CAT_prepare.existing.log'.format(date)
        with open(log_file, 'w') as outf1:
            pass
        
    message = '# CAT v{0}.'.format(about.__version__)
    shared.give_user_feedback(message, log_file, quiet, show_time=False)
    
    message = ('\n'
               'CAT prepare is running, constructing only parts of the '
               'database that are missing. Rawr!\n\n'
               'WARNING: CAT prepare at this point does not check whether the '
               'existing files are OK or corrupted, only if they are there.\n'
               'WARNING: note that the database and taxonomy files should be '
               'downloaded preferably at the same date.\n'
               'WARNING: preparing the database files may take a couple of '
               'hours.\n\n'
               'Supplied command: {0}\n\n'
               'Supplied command: {0}\n\n'
               'Taxonomy folder: {1}/\n'
               'Database folder: {2}/\n'
               'Log file: {3}\n\n'
               '-----------------\n'.format(' '.join(sys.argv),
                                            taxonomy_folder,
                                            database_folder,
                                            log_file))
    shared.give_user_feedback(message, log_file, quiet, show_time=False)
    
    message = 'Doing some pre-flight checks first.'
    shared.give_user_feedback(message, log_file, quiet, show_time=False)
    
    # Check DIAMOND path.
    error = check.check_diamond_binaries(path_to_diamond, log_file, quiet)
    if error:
        sys.exit(1)

    # Check taxonomy folder.
    taxonomy_folder_inspect = check.inspect_taxonomy_folder(taxonomy_folder)
    if taxonomy_folder_inspect == [None]:
        message = ('Taxonomy folder not found. Directory will be created '
                   'fresh and taxonomy files downloaded to it.')
        shared.give_user_feedback(message, log_file, quiet)
        
        nodes_dmp = None
        names_dmp = None
        prot_accession2taxid_file = None
    else:
        (nodes_dmp,
         names_dmp,
         prot_accession2taxid_file) = taxonomy_folder_inspect

        message = ('Taxonomy folder found.')
        shared.give_user_feedback(message, log_file, quiet)
        
    if ((nodes_dmp is None and names_dmp is not None) or
        (nodes_dmp is not None and names_dmp is None)):
        message = ('ERROR: CAT prepare did not find both nodes.dmp and '
                   'names.dmp in the taxonomy folder. They should be '
                   'downloaded together. Remove {0} and try again.'
                   ''.format([file_ for file_ in (nodes_dmp, names_dmp) if
                              file_ is not None][0]))
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    if nodes_dmp is None and names_dmp is None:
        message = ('Nodes.dmp and names.dmp will be downloaded to taxonomy '
                   'folder.')
        shared.give_user_feedback(message, log_file, quiet)

        step_list.append('download_taxonomy_files')
    else:
        message = 'Nodes.dmp found: {0}.'.format(nodes_dmp)
        shared.give_user_feedback(message, log_file, quiet)

        message = 'Names.dmp found: {0}.'.format(names_dmp)
        shared.give_user_feedback(message, log_file, quiet)
        
    if prot_accession2taxid_file is None:
        message = ('Prot.accession2taxid file will be downloaded to taxonomy '
                   'folder.')
        shared.give_user_feedback(message, log_file, quiet)
        
        prot_accession2taxid_file = ('{0}/{1}.prot.accession2taxid.gz'
                                     ''.format(taxonomy_folder, date))
        step_list.append('download_prot_accession2taxid_file')
    else:
        message = ('Prot.accession2taxid file found: {0}.'
                   ''.format(prot_accession2taxid_file))
        shared.give_user_feedback(message, log_file, quiet)

    # Check database folder.
    database_folder_inspect = check.inspect_database_folder(database_folder)
    if database_folder_inspect == [None]:
        message = ('Database folder not found. Directory will be created '
                   'fresh and necessary database files will be downloaded to '
                   'it / constructed in it.')
        shared.give_user_feedback(message, log_file, quiet)
        
        nr_file = None
        diamond_database = None
        fastaid2LCAtaxid_file = None
        taxids_with_multiple_offspring_file = None
    else:
        (nr_file,
         diamond_database,
         fastaid2LCAtaxid_file,
         taxids_with_multiple_offspring_file) = database_folder_inspect

        message = ('Database folder found.')
        shared.give_user_feedback(message, log_file, quiet)

    tmp = (diamond_database,
           fastaid2LCAtaxid_file,
           taxids_with_multiple_offspring_file)
    if (nr_file is None and
        None in tmp and
        not all([file_ is None for file_ in tmp])):
        message = ('ERROR: Database folder does not contain an nr file, while '
                   'some but not all of the downstream files that depend on '
                   'it are present. In order to prevent strange bugs from '
                   'arising, please remove all files from the database folder '
                   'and try again.')
        shared.give_user_feedback(message, log_file, quiet, error=True)
        
        sys.exit(1)

    if (fastaid2LCAtaxid_file is None and
        taxids_with_multiple_offspring_file is not None):
        message = ('ERROR: file taxids_with_multiple_offspring exists but '
                   'fastaid2LCAtaxid is not found in the database folder '
                   'whilst taxids_with_multiple_offspring depends on it. In '
                   'order to prevent strange bugs from arising, please remove '
                   '{0} and try again.'
                   ''.format(taxids_with_multiple_offspring_file))
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)
        
    whether_to_download_nr = True
    if (nr_file is None and
        diamond_database is not None and
        fastaid2LCAtaxid_file is not None and
        taxids_with_multiple_offspring_file is not None):
        whether_to_download_nr = False
        
    if nr_file is None:
        if whether_to_download_nr:
            message = 'Nr file will be downloaded to database folder.'
            shared.give_user_feedback(message, log_file, quiet)
            
            nr_file = '{0}/{1}.nr.gz'.format(database_folder, date)
            step_list.append('download_nr')
        else:
            pass
    else:
        message = 'Nr file found: {0}.'.format(nr_file)
        shared.give_user_feedback(message, log_file, quiet)

    if diamond_database is None:
        message = ('DIAMOND database will be constructed from the nr file.'
                   ''.format(nr_file))
        shared.give_user_feedback(message, log_file, quiet)
        
        diamond_database_prefix = '{0}/{1}.nr'.format(database_folder, date)
        step_list.append('make_diamond_database')
    else:
        message = 'DIAMOND database found: {0}.'.format(diamond_database)
        shared.give_user_feedback(message, log_file, quiet)
        
        diamond_database_prefix = diamond_database.rsplit('.dmnd', 1)[0]

    if fastaid2LCAtaxid_file is None:
        message = 'File fastaid2LCAtaxid will be created.'
        shared.give_user_feedback(message, log_file, quiet)
        
        fastaid2LCAtaxid_file = ('{0}/{1}.nr.fastaid2LCAtaxid'
                                 ''.format(database_folder, date))
        step_list.append('make_fastaid2LCAtaxid_file')
    else:
        message = ('Fastaid2LCAtaxid found: {0}.'
                   ''.format(fastaid2LCAtaxid_file))
        shared.give_user_feedback(message, log_file, quiet)

    if taxids_with_multiple_offspring_file is None:
        message = 'File taxids_with_multiple_offspring will be created.'
        shared.give_user_feedback(message, log_file, quiet)

        taxids_with_multiple_offspring_file = ('{0}/{1}.nr.taxids_with_'
                                               'multiple_offspring'
                                               ''.format(database_folder,
                                                         date))
        step_list.append('make_taxids_with_multiple_offspring_file')
    else:
        message = ('Taxids_with_multiple_offspring found: {0}.'
                   ''.format(taxids_with_multiple_offspring_file))
        shared.give_user_feedback(message, log_file, quiet)

    if nr_file is None and whether_to_download_nr is False:
        # This is pushed here just for the logic of the user.
        message = ('NOTE: Database folder contains all the necessary files '
                   'except for nr.gz. Since nr.gz is not used by CAT or BAT, '
                   'this is fine.')
        shared.give_user_feedback(message, log_file, quiet)

    if taxonomy_folder_inspect == [None] and database_folder_inspect == [None]:
        message = ('\n-----------------\n\n'
                   'WARNING: no taxonomy or database folder was found. CAT '
                   'prepare will create them fresh. Are you sure you are '
                   'linking to existing folders?')
        shared.give_user_feedback(message, log_file, quiet, show_time=False)
        
    if ('make_fastaid2LCAtaxid_file' in step_list or
        'make_taxids_with_multiple_offspring_file' in step_list):
        # Check memory.
        min_mem = 100
        (total_memory, error) = check.check_memory(min_mem)
        
        if error:
            message = ('ERROR: at least {0}GB of memory is needed for the '
                       'database construction. {1}GB is found on your system. '
                       'You can either try to find a machine with more '
                       'memory, or download preconstructed database files '
                       'from tbb.bio.uu.nl/bastiaan/CAT_prepare/.'
                       ''.format(min_mem, total_memory))
            shared.give_user_feedback(message, log_file, quiet, error=True)
            
            sys.exit(1)
            
    if len(step_list) == 0:
        message = ('All necessary files are found. Existing database does not '
                   'need any more work...')
        shared.give_user_feedback(message, log_file, quiet, show_time=False)
    else:
        message = 'Ready to fly!\n\n-----------------\n'
        shared.give_user_feedback(message, log_file, quiet, show_time=False)
        
    if taxonomy_folder_inspect == [None]:
        os.mkdir(taxonomy_folder)
        message = 'Taxonomy folder {0} is created.'.format(taxonomy_folder)
        shared.give_user_feedback(message, log_file, quiet)

    if database_folder_inspect == [None]:
        os.mkdir(database_folder)
        message = 'Database folder {0} is created.'.format(database_folder)
        shared.give_user_feedback(message, log_file, quiet)
        
    prepare(step_list,
            taxonomy_folder,
            database_folder,
            date,
            prot_accession2taxid_file,
            nr_file,
            path_to_diamond,
            diamond_database_prefix,
            nproc,
            fastaid2LCAtaxid_file,
            taxids_with_multiple_offspring_file,
            log_file,
            quiet)
    
def run():
    (args, date) = parse_arguments()
    
    if args.fresh:
        run_fresh(args, date)
    else:
        run_existing(args, date)


if __name__ == '__main__':
    sys.exit('Please run \'CAT prepare\' to construct a CAT/BAT database.')
