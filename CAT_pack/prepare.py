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
    date = datetime.datetime.now().strftime('%Y-%m-%d')
    
    parser = argparse.ArgumentParser(
            prog='CAT prepare',
            description='Download and construct CAT/BAT database files.',
            usage=('CAT prepare (--fresh | --existing) '
                '[options] [-h / --help]'),
            add_help=False)
    
    required_choice = parser.add_argument_group('Required choice')

    group = required_choice.add_mutually_exclusive_group(required=True)

    group.add_argument(
            '--fresh',
            dest='fresh',
            action='store_true',
            help='Start with a fresh database.')
    group.add_argument(
            '--existing',
            dest='fresh',
            action='store_false',
            help=('Start with an existing database. CAT will search the '
                'supplied database and taxonomy folders and only construct '
                'files that do not exist yet.'))
    
    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument(
            '-d',
            '--database_folder',
            dest='database_folder',
            metavar='',
            required=False,
            type=str,
            action=shared.PathAction,
            default='./{0}_CAT_database'.format(date),
            help=('Name of folder to which database files will be written '
                '(default: {date}_CAT_database).'))
    optional.add_argument(
            '-t',
            '--taxonomy_folder',
            dest='taxonomy_folder',
            metavar='',
            required=False,
            type=str,
            action=shared.PathAction,
            default='./{0}_taxonomy'.format(date),
            help=('Name of folder to which taxonomy files will be downloaded '
                '(default: {date}_taxonomy).'))
    optional.add_argument(
            '--path_to_diamond',
            dest='path_to_diamond',
            metavar='',
            required=False,
            type=str,
            action=shared.PathAction,
            default='diamond',
            help=('Path to DIAMOND binaries. Supply if CAT prepare can not '
                'find DIAMOND.'))
    optional.add_argument(
            '-q',
            '--quiet',
            dest='quiet',
            required=False,
            action='store_true',
            help='Suppress verbosity.')
    optional.add_argument(
            '--verbose',
            dest='verbose',
            required=False,
            action='store_true',
            help='Increase verbostity.')
    optional.add_argument(
            '--no_log',
            dest='no_log',
            required=False,
            action='store_true',
            help='Suppress log file.')
    optional.add_argument(
            '-h',
            '--help',
            action='help',
            help='Show this help message and exit.')
    
    specific = parser.add_argument_group('DIAMOND specific optional arguments')

    specific.add_argument(
            '-n',
            '--nproc',
            dest='nproc',
            metavar='',
            required=False,
            type=int,
            default=multiprocessing.cpu_count(),
            help=('Number of cores to deploy by DIAMOND makedb '
                '(default: maximum).'))
    
    (args, extra_args) = parser.parse_known_args()

    extra_args = [arg for (i, arg) in enumerate(extra_args) if
                  (i, arg) != (0, 'prepare')]
    if len(extra_args) > 0:
        sys.exit('error: too much arguments supplied:\n{0}'.format(
            '\n'.join(extra_args)))

    # Add extra arguments.
    setattr(args, 'date', date)
    setattr(args, 'min_mem', 200)
    shared.expand_arguments(args)

    return (args)


def memory_bottleneck(args):
    (total_memory, error) = check.check_memory(args.min_mem)
    if error:
        message = (
                'at least {0}GB of memory is needed for the database '
                'construction. {1}GB is found on your system. You can try '
                'to find a machine with more memory, or download '
                'preconstructed database files from '
                'tbb.bio.uu.nl/bastiaan/CAT_prepare/.'.format(
                    args.min_mem, total_memory))
        shared.give_user_feedback(message, args.log_file, args.quiet,
                error=True)

        sys.exit(1)

    return


def download_taxonomy_files(taxonomy_folder, date, log_file, quiet):
    url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/'
    message = ('Downloading and extracting taxonomy files from {0} to '
            'taxonomy folder.'.format(url))
    shared.give_user_feedback(message, log_file, quiet)

    url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
    tmp_taxonomy_file = '{0}{1}.taxdump.tar.gz'.format(taxonomy_folder, date)
    try:
        urllib.request.urlretrieve(url, tmp_taxonomy_file)
    except:
        message = 'download of {0} failed.'.format(url)
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    url = '{0}.md5'.format(url)
    md5_file = '{0}{1}.taxdump.tar.gz.md5'.format(taxonomy_folder, date)
    try:
        urllib.request.urlretrieve(url, md5_file)
    except:
        message = 'download of {0} failed.'.format(url)
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    message = 'Download complete.'
    shared.give_user_feedback(message, log_file, quiet)

    check.check_md5_gz(tmp_taxonomy_file, md5_file, log_file, quiet)

    try:
        with tarfile.open(tmp_taxonomy_file) as tar:
            tar.extractall(taxonomy_folder)
    except:
        message = 'something went wrong while extracting the taxonomy files.'
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    message = 'Extracting complete.'
    shared.give_user_feedback(message, log_file, quiet)

    return
    
    
def download_prot_accession2taxid_file(
        prot_accession2taxid_file, date, log_file, quiet):
    url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/'
    message = 'Downloading mapping file from {0} to taxonomy folder.'.format(
            url)
    shared.give_user_feedback(message, log_file, quiet)

    url = ('ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/'
            'prot.accession2taxid.FULL.gz')
    try:
        urllib.request.urlretrieve(url, prot_accession2taxid_file)
    except:
        message = 'download of {0} failed.'.format(url)
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    url = '{0}.md5'.format(url)
    md5_file = '{0}.md5'.format(prot_accession2taxid_file)
    try:
        urllib.request.urlretrieve(url, md5_file)
    except:
        message = 'download of {0} failed.'.format(url)
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    message = 'Download complete.'
    shared.give_user_feedback(message, log_file, quiet)

    check.check_md5_gz(prot_accession2taxid_file, md5_file, log_file, quiet)

    return


def download_nr(nr_file, log_file, quiet):
    url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/'
    message = 'Downloading nr database from {0} to database folder.'.format(
            url)
    shared.give_user_feedback(message, log_file, quiet)

    url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz'
    try:
        urllib.request.urlretrieve(url, nr_file)
    except:
        message = 'download of {0} failed.'.format(url)
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    url = '{0}.md5'.format(url)
    md5_file = '{0}.md5'.format(nr_file)
    try:
        urllib.request.urlretrieve(url, md5_file)
    except:
        message = 'download of {0} failed.'.format(url)
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    message = 'Download complete.'
    shared.give_user_feedback(message, log_file, quiet)

    check.check_md5_gz(nr_file, md5_file, log_file, quiet)

    return


def make_diamond_database(
        path_to_diamond,
        nr_file,
        diamond_database_prefix,
        nproc,
        log_file,
        quiet,
        verbose):
    message = ('Constructing DIAMOND database {0}.dmnd from {1} using {2} '
            'cores.'.format(diamond_database_prefix, nr_file, nproc))
    shared.give_user_feedback(message, log_file, quiet)

    command = [
            path_to_diamond, 'makedb',
            '--in', nr_file,
            '-d', diamond_database_prefix,
            '-p', str(nproc)]

    if not verbose:
        command += ['--quiet']

    try:
        subprocess.check_call(command)
    except:
        message = 'DIAMOND database could not be created.'
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)
        
    message = 'DIAMOND database constructed.'
    shared.give_user_feedback(message, log_file, quiet)

    return


def import_headers_nr(nr_file, log_file, quiet):
    message = 'Loading file {0}.'.format(nr_file)
    shared.give_user_feedback(message, log_file, quiet)

    fastaid2prot_accessions = {}
    prot_accessions_whitelist = set()

    with gzip.open(nr_file, 'rb') as f1:
        for line in f1:
            line = line.decode('utf-8')

            if not line.startswith('>'):
                continue

            line = line.lstrip('>').split('\x01')

            prot_accessions = [i.split(' ')[0] for i in line]
            fastaid = prot_accessions[0]

            fastaid2prot_accessions[fastaid] = prot_accessions
            prot_accessions_whitelist.update(prot_accessions)

    return (fastaid2prot_accessions, prot_accessions_whitelist)


def import_prot_accession2taxid(
        prot_accession2taxid_file, prot_accessions_whitelist, log_file, quiet):
    message = 'Loading file {0}.'.format(prot_accession2taxid_file)
    shared.give_user_feedback(message, log_file, quiet)
    
    prot_accession2taxid = {}

    with gzip.open(prot_accession2taxid_file, 'rb') as f1:
        for n, line in enumerate(f1):
            line = line.decode('utf-8')

            line = line.rstrip().split('\t')

            if n == 0:
                index_1 = line.index('accession.version')
                index_2 = line.index('taxid')

                continue

            prot_accession = line[index_1]

            if prot_accession in prot_accessions_whitelist:
                prot_accession2taxid[prot_accession] = line[index_2]

    return prot_accession2taxid


def make_fastaid2LCAtaxid_file(
        nodes_dmp,
        fastaid2LCAtaxid_file,
        nr_file,
        prot_accession2taxid_file,
        taxid2parent,
        log_file,
        quiet):
    (fastaid2prot_accessions,
            prot_accessions_whitelist) = import_headers_nr(
                    nr_file, log_file, quiet)
    prot_accession2taxid = import_prot_accession2taxid(
            prot_accession2taxid_file, prot_accessions_whitelist,
            log_file, quiet)

    message = 'Finding LCA of all protein accession numbers in fasta headers.'
    shared.give_user_feedback(message, log_file, quiet)
    
    no_taxid = 0
    corrected = 0
    total = 0
    with open(fastaid2LCAtaxid_file, 'w') as outf1:
        for fastaid, prot_accessions in fastaid2prot_accessions.items():
            list_of_lineages = []
            for prot_accession in prot_accessions:
                try:
                    taxid = prot_accession2taxid[prot_accession]
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
                # missing in nodes.dmp. NOTE that these entries are thus not
                # present in the output file.
                no_taxid += 1

                continue

            LCAtaxid = tax.find_LCA(list_of_lineages)

            outf1.write('{0}\t{1}\n'.format(fastaid, LCAtaxid))

            if (fastaid not in prot_accession2taxid or
                    LCAtaxid != prot_accession2taxid[fastaid]):
                # If the fastaid cannot be found in prot.accession2taxid, but
                # a taxid is given to the fastaid based on secondary accession
                # numbers, or if the taxid of the header is different from the
                # LCA taxid, it is counted as corrected.
                corrected += 1

    message = ('Done! File {0} is created. '
            '{1:,d} of {2:,d} headers ({3:.1f}%) corrected. '
            '{4:,d} headers ({5:.1f}%) do not have a taxid assigned.'.format(
                fastaid2LCAtaxid_file,
                corrected,
                total,
                corrected / total * 100,
                no_taxid,
                no_taxid / total * 100))
    shared.give_user_feedback(message, log_file, quiet)

    return


def find_offspring(
        nodes_dmp, fastaid2LCAtaxid_file, taxid2parent, log_file, quiet):
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


def write_taxids_with_multiple_offspring_file(
        taxids_with_multiple_offspring_file, taxid2offspring, log_file, quiet):
    message = 'Writing {0}.'.format(taxids_with_multiple_offspring_file)
    shared.give_user_feedback(message, log_file, quiet)

    with open(taxids_with_multiple_offspring_file, 'w') as outf1:
        for taxid in taxid2offspring:
            if len(taxid2offspring[taxid]) >= 2:
                outf1.write('{0}\n'.format(taxid))

    return
                        
                        
def prepare(step_list, args):
    shared.print_variables(args, step_list)

    if not os.path.isdir(args.taxonomy_folder):
        os.mkdir(args.taxonomy_folder)
        message = 'Taxonomy folder {0} is created.'.format(
                args.taxonomy_folder)
        shared.give_user_feedback(message, args.log_file, args.quiet)

    if not os.path.isdir(args.database_folder):
        os.mkdir(args.database_folder)
        message = 'Database folder {0} is created.'.format(
                args.database_folder)
        shared.give_user_feedback(message, args.log_file, args.quiet)

    if 'download_taxonomy_files' in step_list:
        download_taxonomy_files(
                args.taxonomy_folder, args.date, args.log_file, args.quiet)

        setattr(args, 'nodes_dmp', '{0}nodes.dmp'.format(args.taxonomy_folder))

    if 'download_prot_accession2taxid_file' in step_list:
        setattr(args,
                'prot_accession2taxid_file',
                '{0}{1}.prot.accession2taxid.FULL.gz'.format(
                    args.taxonomy_folder, args.date))

        download_prot_accession2taxid_file(
                args.prot_accession2taxid_file,
                args.date,
                args.log_file,
                args.quiet)

    if 'download_nr' in step_list:
        setattr(args,
                'nr_file',
                '{0}{1}.nr.gz'.format(args.database_folder, args.date))

        download_nr(args.nr_file, args.log_file, args.quiet)

    if 'make_diamond_database' in step_list:
        setattr(args,
                'diamond_database_prefix',
                '{0}{1}.nr'.format(args.database_folder, args.date))

        make_diamond_database(
                args.path_to_diamond,
                args.nr_file,
                args.diamond_database_prefix,
                args.nproc,
                args.log_file,
                args.quiet,
                args.verbose)

    if ('make_fastaid2LCAtaxid_file' in step_list
            or 'make_taxids_with_multiple_offspring_file' in step_list):
        taxid2parent, taxid2rank = tax.import_nodes(
                args.nodes_dmp, args.log_file, args.quiet)

    if 'make_fastaid2LCAtaxid_file' in step_list:
        setattr(args,
                'fastaid2LCAtaxid_file',
                '{0}{1}.nr.fastaid2LCAtaxid'.format(
                    args.database_folder, args.date))

        make_fastaid2LCAtaxid_file(
                args.nodes_dmp,
                args.fastaid2LCAtaxid_file,
                args.nr_file,
                args.prot_accession2taxid_file,
                taxid2parent,
                args.log_file,
                args.quiet)

    if 'make_taxids_with_multiple_offspring_file' in step_list:
        setattr(args,
                'taxids_with_multiple_offspring_file',
                '{0}{1}.nr.taxids_with_multiple_offspring'.format(
                    args.database_folder, args.date))

        taxid2offspring = find_offspring(
                args.nodes_dmp,
                args.fastaid2LCAtaxid_file,
                taxid2parent,
                args.log_file,
                args.quiet)
        write_taxids_with_multiple_offspring_file(
                args.taxids_with_multiple_offspring_file,
                taxid2offspring,
                args.log_file,
                args.quiet)

    message = ('\n-----------------\n\n'
            '{0} CAT prepare is done!'.format(shared.timestamp()))
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)

    if args.nr_file:
        message = 'You may remove {0} now.'.format(args.nr_file)
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=False)

    message = (
            '\nSupply the following arguments to CAT or BAT if you want to '
            'use this database:\n'
            '-d / --database_folder {0}\n'
            '-t / --taxonomy_folder {1}'.format(
                args.database_folder, args.taxonomy_folder))
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)

    return
    
    
def run_fresh(args):
    message = '# CAT v{0}.'.format(about.__version__)
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)
    
    message = (
            '\n'
            'CAT prepare is running, constructing a fresh database.\n'
            'Rawr!\n\n'
            'WARNING: preparing the database files may take a couple of hours.'
            '\n\n'
            'Supplied command: {0}\n\n'
            'Taxonomy folder: {1}\n'
            'Database folder: {2}\n'
            'Log file: {3}\n\n'
            '-----------------\n'.format(
                ' '.join(sys.argv),
                args.taxonomy_folder,
                args.database_folder,
                args.log_file))
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)
    
    # Check diamond path.
    error = check.check_diamond_binaries(
            args.path_to_diamond, args.log_file, args.quiet)
    if error:
        sys.exit(1)

    if os.path.isdir(args.taxonomy_folder):
        if args.nodes_dmp or args.names_dmp or args.prot_accession2taxid_file:
            message = (
                    'taxonomy folder {0} exists already and contains taxonomy '
                    'files. Supply a novel or empty folder if you want '
                    'to start fresh, or run CAT prepare --existing.'.format(
                        args.taxonomy_folder))
            shared.give_user_feedback(message, args.log_file, args.quiet,
                    error=True)

            sys.exit(1)

        message = ('Taxonomy folder exists already. Taxonomy files will be '
                'downloaded to it.')
        shared.give_user_feedback(message, args.log_file, args.quiet)

    if os.path.isdir(args.database_folder):
        if (args.nr_file or
                args.diamond_database or
                args.fastaid2LCAtaxid_file or
                args.taxids_with_multiple_offspring_file):
            message = (
                    'database folder {0} exists already and contains database '
                    'files. Supply a novel or empty folder if you want to '
                    'start fresh.'.format(args.database_folder))
            shared.give_user_feedback(message, args.log_file, args.quiet,
                    error=True)

            sys.exit(1)

        message = ('Database folder exists already. Database file will be '
                'downloaded to it / constructed in it.')
        shared.give_user_feedback(message, args.log_file, args.quiet)
        
    # Check memory.
    memory_bottleneck(args)

    step_list = ['download_taxonomy_files',
                 'download_prot_accession2taxid_file',
                 'download_nr',
                 'make_diamond_database',
                 'make_fastaid2LCAtaxid_file',
                 'make_taxids_with_multiple_offspring_file']
    
    prepare(step_list, args)

    return
    
    
def run_existing(args):
    step_list = []

    message = '# CAT v{0}.'.format(about.__version__)
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)
    
    message = (
            '\n'
            'CAT prepare is running, constructing only parts of the database '
            'that are missing. Rawr!\n\n'
            'WARNING: CAT prepare does not check whether the existing files '
            'are OK or corrupted, only if they are there.\n'
            'WARNING: note that the database and taxonomy files should be '
            'downloaded preferably at the same date.\n'
            'WARNING: preparing the database files may take a couple of hours.'
            '\n\n'
            'Supplied command: {0}\n\n'
            'Taxonomy folder: {1}\n'
            'Database folder: {2}\n'
            'Log file: {3}\n\n'
            '-----------------\n'.format(
                ' '.join(sys.argv),
                args.taxonomy_folder,
                args.database_folder,
                args.log_file))
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)
    
    message = 'Doing some pre-flight checks first.'
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)
    
    # Check DIAMOND path.
    error = check.check_diamond_binaries(
            args.path_to_diamond, args.log_file, args.quiet)
    if error:
        sys.exit(1)

    # Check taxonomy folder.
    if not os.path.isdir(args.taxonomy_folder):
        message = ('Taxonomy folder not found. Directory will be created '
                'fresh and taxonomy files downloaded to it.')
        shared.give_user_feedback(message, args.log_file, args.quiet)
    else:
        message = ('Taxonomy folder found.')
        shared.give_user_feedback(message, args.log_file, args.quiet)

    if ((not args.nodes_dmp and args.names_dmp) or
            (args.nodes_dmp and not args.names_dmp)):
        message = (
                'CAT prepare did not find both nodes.dmp and names.dmp in the '
                'taxonomy folder. They should be downloaded together. Remove '
                '{0} and try again.'.format(
                    [file_ for file_ in (args.nodes_dmp, args.names_dmp) if
                        file_][0]))
        shared.give_user_feedback(message, args.log_file, args.quiet,
                error=True)

        sys.exit(1)

    if not args.nodes_dmp and not args.names_dmp:
        message = ('Nodes.dmp and names.dmp will be downloaded to taxonomy '
                'folder.')
        shared.give_user_feedback(message, args.log_file, args.quiet)

        step_list.append('download_taxonomy_files')
    else:
        message = 'Nodes.dmp found: {0}.'.format(args.nodes_dmp)
        shared.give_user_feedback(message, args.log_file, args.quiet)

        message = 'Names.dmp found: {0}.'.format(args.names_dmp)
        shared.give_user_feedback(message, args.log_file, args.quiet)

    if not args.prot_accession2taxid_file:
        # NOTE that the file will only be downloaded if a new
        # fastaid2LCAtaxid_file needs to be constructed.
        message = 'Prot.accession2taxid file not found in taxonomy folder.'
        shared.give_user_feedback(message, args.log_file, args.quiet)
    else:
        message = 'Prot.accession2taxid file found: {0}.'.format(
                args.prot_accession2taxid_file)
        shared.give_user_feedback(message, args.log_file, args.quiet)

    # Check database folder.
    if not os.path.isdir(args.database_folder):
        message = (
                'Database folder not found. Directory will be created fresh '
                'and necessary database files will be downloaded to '
                'it / constructed in it.')
        shared.give_user_feedback(message, args.log_file, args.quiet)
    else:
        message = ('Database folder found.')
        shared.give_user_feedback(message, args.log_file, args.quiet)

    tmp = (args.diamond_database,
            args.fastaid2LCAtaxid_file,
            args.taxids_with_multiple_offspring_file)

    if (not args.nr_file and
            None in tmp and
            not all([file_ is None for file_ in tmp])):
        message = (
                'database folder does not contain an nr file, while some but '
                'not all of the downstream files that depend on it are '
                'present. In order to prevent strange bugs from arising, '
                'remove all files from the database folder and try again.')
        shared.give_user_feedback(message, args.log_file, args.quiet,
                error=True)
        
        sys.exit(1)

    if (not args.fastaid2LCAtaxid_file and
            args.taxids_with_multiple_offspring_file):
        message = (
                'file taxids_with_multiple_offspring exists but '
                'fastaid2LCAtaxid is not found in the database folder whilst '
                'taxids_with_multiple_offspring depends on it. In order to '
                'prevent strange bugs from arising, remove {0} and try again.'
                ''.format(args.taxids_with_multiple_offspring_file))
        shared.give_user_feedback(message, args.log_file, args.quiet,
                error=True)

        sys.exit(1)
        
    whether_to_download_nr = True
    if (not args.nr_file and
        args.diamond_database and
        args.fastaid2LCAtaxid_file and
        args.taxids_with_multiple_offspring_file):
        whether_to_download_nr = False
        
    if not args.nr_file:
        if whether_to_download_nr:
            message = 'Nr file will be downloaded to database folder.'
            shared.give_user_feedback(message, args.log_file, args.quiet)

            step_list.append('download_nr')
        else:
            pass
    else:
        message = 'Nr file found: {0}.'.format(args.nr_file)
        shared.give_user_feedback(message, args.log_file, args.quiet)

    if not args.diamond_database:
        message = 'DIAMOND database will be constructed from the nr file.'
        shared.give_user_feedback(message, args.log_file, args.quiet)

        step_list.append('make_diamond_database')
    else:
        message = 'DIAMOND database found: {0}.'.format(args.diamond_database)
        shared.give_user_feedback(message, args.log_file, args.quiet)

    if not args.fastaid2LCAtaxid_file:
        if not args.prot_accession2taxid_file:
            message = ('Prot.accession2taxid file will be downloaded to '
                    'taxonomy folder.')
            shared.give_user_feedback(message, args.log_file, args.quiet)

            step_list.append('download_prot_accession2taxid_file')

        message = 'File fastaid2LCAtaxid will be created.'
        shared.give_user_feedback(message, args.log_file, args.quiet)

        step_list.append('make_fastaid2LCAtaxid_file')
    else:
        message = ('Fastaid2LCAtaxid found: {0}.'.format(
            args.fastaid2LCAtaxid_file))
        shared.give_user_feedback(message, args.log_file, args.quiet)

        if not args.prot_accession2taxid_file:
            message = 'Prot.accession2taxid file will not be needed.'
            shared.give_user_feedback(message, args.log_file, args.quiet)

    if not args.taxids_with_multiple_offspring_file:
        message = 'File taxids_with_multiple_offspring will be created.'
        shared.give_user_feedback(message, args.log_file, args.quiet)

        step_list.append('make_taxids_with_multiple_offspring_file')
    else:
        message = 'Taxids_with_multiple_offspring found: {0}'.format(
            args.taxids_with_multiple_offspring_file)
        shared.give_user_feedback(message, args.log_file, args.quiet)

    if not args.nr_file and whether_to_download_nr is False:
        # This is pushed here just for the logic of the user.
        message = (
                'NOTE: Database folder contains all the necessary files '
                'except for nr.gz. Since nr.gz is not used by CAT or BAT, '
                'this is fine.')
        shared.give_user_feedback(message, args.log_file, args.quiet)

    if (not os.path.isdir(args.taxonomy_folder) and
            not os.path.isdir(args.database_folder)):
        message = (
                '\n-----------------\n\n'
                'WARNING: no taxonomy or database folder was found. CAT '
                'prepare will create them fresh. Are you sure you are linking '
                'to existing folders?')
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=False)

    if 'make_fastaid2LCAtaxid_file' in step_list:
        # Check memory.
        memory_bottleneck(args)

    if len(step_list) == 0:
        message = ('All necessary files are found. Existing database does not '
                'need any more work...')
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=False)

        sys.exit(0)
    else:
        message = 'Ready to fly!\n\n-----------------\n'
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=False)

    prepare(step_list, args)

    return

    
def run():
    args = parse_arguments()
    
    if args.fresh:
        run_fresh(args)
    else:
        run_existing(args)

    return


if __name__ == '__main__':
    sys.exit('Run \'CAT prepare\' to construct a CAT/BAT database.')
