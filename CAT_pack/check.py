#!/usr/bin/env/ python3

import datetime
import os
import subprocess
import sys

import shared


def convert_arguments(args):
    if 'named_input_file' in args:
        # A call from summarise.
        return (args.named_input_file,
                args.output_file,
                args.contigs_fasta,
                args.force,
                args.quiet)
        
    taxonomy_folder = os.path.expanduser(args.taxonomy_folder.rstrip('/'))

    if 'only_official' in args:
        # A call from add_names.
        return (args.input_file,
                args.output_file,
                taxonomy_folder,
                args.only_official,
                args.force,
                args.quiet)
        
    database_folder = os.path.expanduser(args.database_folder.rstrip('/'))

    if 'fresh' in args:
        # A call from prepare.
        return (database_folder,
                taxonomy_folder,
                args.path_to_diamond,
                args.quiet,
                args.no_log,
                args.nproc)
    
    one_minus_r = (100 - args.r) / 100

    if args.tmpdir is None:
        if '/' in args.out_prefix:
            tmpdir = args.out_prefix.rsplit('/', 1)[0]
        else:
            tmpdir = './'
    else:
        tmpdir = args.tmpdir
        
    if 'bin_suffix' in args:
        # A call from bins.
        bin_folder = os.path.expanduser(args.bin_folder.rstrip('/'))
        bin_suffix = '.{0}'.format(args.bin_suffix.lstrip('.'))
        
        return (bin_folder,
                database_folder,
                taxonomy_folder,
                bin_suffix,
                one_minus_r,
                args.f,
                args.out_prefix,
                args.predicted_proteins_fasta,
                args.diamond_file,
                args.path_to_prodigal,
                args.path_to_diamond,
                args.force,
                args.quiet,
                args.no_log,
                args.nproc,
                args.sensitive,
                args.block_size,
                args.index_chunks,
                tmpdir)
    else:
        # A call from contigs.
        return (args.contigs_fasta,
                database_folder,
                taxonomy_folder,
                one_minus_r,
                args.f,
                args.out_prefix,
                args.predicted_proteins_fasta,
                args.diamond_file,
                args.path_to_prodigal,
                args.path_to_diamond,
                args.force,
                args.quiet,
                args.no_log,
                args.nproc,
                args.sensitive,
                args.block_size,
                args.index_chunks,
                tmpdir)


def check_memory(Gb):
    error = False
    
    if sys.platform == 'linux' or sys.platform == 'linux2':
        # It's a Linux!
        meminfo_file = '/proc/meminfo'
        with open(meminfo_file, 'r') as f:
            for line in f:
                if line.startswith('MemTotal:'):
                    mem = int(line.split(' ')[-2])

                    # Mem is given in Kb, convert to Gb.
                    total_memory = mem / 2 ** 20
    elif sys.platform == 'darwin':
        # It's a Mac!
        meminfo = subprocess.check_output(['sysctl', 'hw.memsize'])
        mem = int(meminfo.decode('utf-8').rstrip().split(' ')[-1])
        
        # Mem is given in b, convert to Gb.
        total_memory = mem / 2 ** 30
        
    if total_memory < Gb:
        error = True
        
    return ('{0:.1f}'.format(total_memory), error)


def check_out_prefix(out_prefix, log_file, quiet):
    error = False

    if '/' in out_prefix:
        if out_prefix.endswith('/'):
            message = ('ERROR: prefix for output files ({0}) appears to be a '
                       'directory.'.format(out_prefix))
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

        directory = out_prefix.rsplit('/', 1)[0]
        if not os.path.isdir(directory):
            message = ('ERROR: can not find output directory {0} to which '
                       'output files should be written.'.format(directory))
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

    return error


def check_diamond_binaries(path_to_diamond, log_file, quiet):
    error = False

    try:
        p = subprocess.Popen([path_to_diamond, '--version'],
                             stdout=subprocess.PIPE)
        c = p.communicate()
        output = c[0].decode().rstrip()

        message = 'DIAMOND found: {0}.'.format(output)
        shared.give_user_feedback(message, log_file, quiet)
    except OSError:
        message = ('ERROR: can not find DIAMOND. Please check whether it is '
                   'installed or path to the binaries is provided.')
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_prodigal_binaries(path_to_prodigal, log_file, quiet):
    error = False

    try:
        p = subprocess.Popen([path_to_prodigal, '-v'], stderr=subprocess.PIPE)
        c = p.communicate()
        output = c[1].decode().rstrip().lstrip()

        message = 'Prodigal found: {0}.'.format(output)
        shared.give_user_feedback(message, log_file, quiet)
    except OSError:
        message = ('ERROR: can not find Prodigal. Please check whether it is '
                   'installed or path to the binaries is provided.')
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_bin_folder(bin_folder, bin_suffix, log_file, quiet):
    error = False

    if not os.path.isdir(bin_folder):
        message = ('ERROR: can not find the bin folder.')
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

        return error
    
    tmp = []
    for file in os.listdir(bin_folder):
        if file.endswith(bin_suffix):
            tmp.append(file)

    if len(tmp) == 0:
        message = ('ERROR: no bins found with suffix {0} in bin folder. You '
                   'can set the suffix with the [-s / --bin_suffix] argument.'
                   ''.format(bin_suffix))
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def inspect_taxonomy_folder(taxonomy_folder):
    if not os.path.isdir(taxonomy_folder):
        return [None]

    nodes_dmp = None 
    names_dmp = None
    prot_accession2taxid_file = None

    for file in os.listdir(taxonomy_folder):
        if file == 'nodes.dmp':
            nodes_dmp = '{0}/{1}'.format(taxonomy_folder, file)
        elif file == 'names.dmp':
            names_dmp = '{0}/{1}'.format(taxonomy_folder, file)
        elif file.endswith('prot.accession2taxid.gz'):
            prot_accession2taxid_file = '{0}/{1}'.format(taxonomy_folder, file)

    return (nodes_dmp, names_dmp, prot_accession2taxid_file)


def inspect_database_folder(database_folder):
    if not os.path.isdir(database_folder):
        return [None]

    nr_file = None
    diamond_database = None
    fastaid2LCAtaxid_file = None
    taxids_with_multiple_offspring_file = None

    for file in os.listdir(database_folder):
        if file.endswith('nr.gz'):
            nr_file = '{0}/{1}'.format(database_folder, file)
        elif file.endswith('.dmnd'):
            diamond_database = '{0}/{1}'.format(database_folder, file)
        elif file.endswith('fastaid2LCAtaxid'):
            fastaid2LCAtaxid_file = '{0}/{1}'.format(database_folder, file)
        elif file.endswith('taxids_with_multiple_offspring'):
            taxids_with_multiple_offspring_file = ('{0}/{1}'
                                                   ''.format(database_folder,
                                                             file))
                                                   
    return (nr_file,
            diamond_database,
            fastaid2LCAtaxid_file,
            taxids_with_multiple_offspring_file)
    
    
def check_folders_for_run(taxonomy_folder,
                          database_folder,
                          step_list,
                          log_file, quiet):
    error = False

    taxonomy_folder_inspect = inspect_taxonomy_folder(taxonomy_folder)
    if taxonomy_folder_inspect == [None]:
        message = 'ERROR: can not find the taxonomy folder.'
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True
    else:
        (nodes_dmp,
         names_dmp,
         prot_accession2taxid_file) = taxonomy_folder_inspect

        if nodes_dmp is None or names_dmp is None:
            message = ('ERROR: nodes.dmp and / or names.dmp not found in the '
                       'taxonomy folder.')
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

    database_folder_inspect = inspect_database_folder(database_folder)
    if database_folder_inspect == [None]:
        message = 'ERROR: can not find the database folder.'
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True
    else:
        (nr_file,
         diamond_database,
         fastaid2LCAtaxid_file,
         taxids_with_multiple_offspring_file) = database_folder_inspect

        if diamond_database is None and 'run_diamond' in step_list:
            message = 'ERROR: DIAMOND database not found in database folder.'
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

        if fastaid2LCAtaxid_file is None:
            message = ('ERROR: file fastaid2LCAtaxid is not found in database '
                       'folder.' )
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

        if taxids_with_multiple_offspring_file is None:
            message = ('ERROR: file taxids_with_multiple_offspring not found '
                       'in database folder.')
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

    return error


def check_output_file(output_file, log_file, quiet):
    error = False

    if os.path.isfile(output_file):
        message = ('ERROR: output file {0} already exists. You can choose to '
                   'overwrite existing files with the [--force] argument.'
                   ''.format(output_file))
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_input_file(input_file, log_file, quiet):
    error = False

    if not os.path.isfile(input_file):
        message = 'ERROR: input file {0} does not exist.'.format(input_file)
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_whether_file_is_fasta(file):
    is_fasta = False

    if not os.path.isfile(file):
        return is_fasta
    
    with open(file, 'r') as f1:
        for line in f1:
            if line.startswith('>'):
                is_fasta = True

            break

    return is_fasta


def check_whether_ORFs_are_based_on_contigs(contig_names,
                                            contig2ORFs,
                                            log_file,
                                            quiet):
    for contig in contig2ORFs:
        if contig not in contig_names:
            message = ('ERROR: found a protein in the predicted proteins '
                       'fasta file that can not be traced back to one of the '
                       'contigs in the contigs fasta file: {0}. Proteins '
                       'should be named contig_name_#.'
                       ''.format(contig2ORFs[contig][0]))
            shared.give_user_feedback(message, log_file, quiet, error=True)

            sys.exit(1)
            
            
if __name__ == '__main__':
    sys.exit('Please run \'CAT\' to run CAT or BAT.')
