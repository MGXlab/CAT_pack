#!/usr/bin/env python3

import argparse
import datetime
import decimal
import gzip
import os
import subprocess
import sys

import check


class PathAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        path = os.path.expanduser(values.rstrip('/'))

        if not path.startswith('/') and not path.startswith('.'):
            path = './{0}'.format(path)

        if os.path.isdir(path):
            path = '{0}/'.format(path)

        setattr(namespace, self.dest, path)


class DecimalAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, decimal.Decimal(values))


class SuffixAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        bin_suffix = '.{0}'.format(values.lstrip('.'))

        setattr(namespace, self.dest, bin_suffix)


def expand_arguments(args):
    if 'r' in args:
        setattr(args, 'one_minus_r', (100 - args.r) / 100)

    if 'out_prefix' in args:
        if not args.tmpdir:
            tmpdir = '{0}/'.format(args.out_prefix.rsplit('/', 1)[0])

            setattr(args, 'tmpdir', tmpdir)

    if 'no_log' in args and not args.no_log:
        if 'fresh' in args and args.fresh:
            log_file = './{0}.CAT_prepare.fresh.log'.format(args.date)
        elif 'fresh' in args and not args.fresh:
            log_file = './{0}.CAT_prepare.existing.log'.format(args.date)
        else:
            # Check out_prefix as the log file needs to be written to a valid
            # location.
            error = check.check_out_prefix(args.out_prefix, None, args.quiet)
            if error:
                sys.exit(1)

            log_file = '{0}.log'.format(args.out_prefix)

        with open(log_file, 'w') as outf:
            pass
    else:
        log_file = None

    setattr(args, 'log_file', log_file)

    if 'taxonomy_folder' in args:
        setattr(args,
                'taxonomy_folder',
                '{0}/'.format(args.taxonomy_folder.rstrip('/')))

        explore_taxonomy_folder(args)
    if 'database_folder' in args:
        setattr(args,
                'database_folder',
                '{0}/'.format(args.database_folder.rstrip('/')))

        explore_database_folder(args)

    return


def explore_taxonomy_folder(args):
    nodes_dmp = None
    names_dmp = None
    prot_accession2taxid_file = None

    if os.path.isdir(args.taxonomy_folder):
        for file_ in os.listdir(args.taxonomy_folder):
            if file_ == 'nodes.dmp':
                nodes_dmp = '{0}{1}'.format(args.taxonomy_folder, file_)
            elif file_ == 'names.dmp':
                names_dmp = '{0}{1}'.format(args.taxonomy_folder, file_)
            elif file_.endswith('prot.accession2taxid.gz'):
                prot_accession2taxid_file = '{0}{1}'.format(
                        args.taxonomy_folder, file_)

    setattr(args, 'nodes_dmp', nodes_dmp)
    setattr(args, 'names_dmp', names_dmp)
    setattr(args, 'prot_accession2taxid_file', prot_accession2taxid_file)

    return


def explore_database_folder(args):
    nr_file = None
    diamond_database = None
    fastaid2LCAtaxid_file = None
    taxids_with_multiple_offspring_file = None

    if os.path.isdir(args.database_folder):
        for file_ in os.listdir(args.database_folder):
            if file_.endswith('nr.gz'):
                nr_file = '{0}{1}'.format(args.database_folder, file_)
            elif file_.endswith('.dmnd'):
                diamond_database = '{0}{1}'.format(
                        args.database_folder, file_)
            elif file_.endswith('fastaid2LCAtaxid'):
                fastaid2LCAtaxid_file = '{0}{1}'.format(
                        args.database_folder, file_)
            elif file_.endswith('taxids_with_multiple_offspring'):
                taxids_with_multiple_offspring_file = ('{0}{1}'
                        ''.format(args.database_folder, file_))

    setattr(args, 'nr_file', nr_file)
    setattr(args, 'diamond_database', diamond_database)
    setattr(args, 'fastaid2LCAtaxid_file', fastaid2LCAtaxid_file)
    setattr(args,
            'taxids_with_multiple_offspring_file',
            taxids_with_multiple_offspring_file)

    return


def print_variables(args, step_list=None):
    if args.verbose:
        arguments = ['{0}: {1}'.format(k, v) for
                k, v in sorted(vars(args).items())]
        message = (
                '\n-----------------\n\n'
                'Full list of arguments:\n'
                '{0}'.format('\n'.join(arguments)))
        give_user_feedback(message, args.log_file, args.quiet, show_time=False)

        if step_list is not None:
            message = '\nStep list: {0}'.format(step_list)
            give_user_feedback(message, args.log_file, args.quiet,
                    show_time=False)

        message = '\n-----------------\n'
        give_user_feedback(message, args.log_file, args.quiet, show_time=False)

    return


def give_user_feedback(
        message, log_file=None, quiet=False, show_time=True, error=False):
    if error:
        message = 'ERROR: {0}'.format(message)

    if show_time:
        time = datetime.datetime.now()

        message = '[{0}] {1}'.format(time, message)

    message = '{0}\n'.format(message)

    if log_file:
        with open(log_file, 'a') as outf1:
            outf1.write(message)

    if not quiet and not error:
        sys.stdout.write(message)

    if not quiet and error:
        sys.stderr.write(message)

    return
        
        
def run_prodigal(
        path_to_prodigal,
        contigs_fasta,
        proteins_fasta,
        proteins_gff,
        log_file,
        quiet):
    message = (
            'Running Prodigal for ORF prediction. Files {0} and {1} will be '
            'generated. Do not forget to cite Prodigal when using CAT or BAT '
            'in your publication!'.format(proteins_fasta, proteins_gff))
    give_user_feedback(message, log_file, quiet)

    try:
        command = [
                path_to_prodigal,
                '-i', contigs_fasta,
                '-a', proteins_fasta,
                '-o', proteins_gff,
                '-p', 'meta',
                '-g', '11',
                '-q',
                '-f', 'gff']
        subprocess.check_call(command)
    except:
        message = 'Prodigal finished abnormally.'
        give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    message = 'ORF prediction done!'
    give_user_feedback(message, log_file, quiet)

    return
    
    
def run_diamond(args):
    if args.sensitive:
        mode = 'sensitive'
    else:
        mode = 'fast'

    if args.compress:
        compression = '1'
    else:
        compression = '0'

    message = (
            'Homology search with DIAMOND is starting. Please be patient. Do '
            'not forget to cite DIAMOND when using CAT or BAT in your '
            'publication!\n'
            '\t\t\t\tquery: {0}\n'
            '\t\t\t\tdatabase: {1}\n'
            '\t\t\t\tmode: {2}\n'
            '\t\t\t\tnumber of cores: {3}\n'
            '\t\t\t\tblock-size (billions of letters): {4}\n'
            '\t\t\t\tindex-chunks: {5}\n'
            '\t\t\t\ttmpdir: {6}\n'
            '\t\t\t\tcompress: {7}\n'
            '\t\t\t\ttop: {8}'.format(
                args.proteins_fasta,
                args.diamond_database,
                mode,
                args.nproc,
                args.block_size,
                args.index_chunks,
                args.tmpdir,
                compression,
                args.top))
    give_user_feedback(message, args.log_file, args.quiet)

    try:
        command = [
                args.path_to_diamond,
                'blastp',
                '-d', args.diamond_database,
                '-q', args.proteins_fasta,
                '--top', str(args.top),
                '--matrix', 'BLOSUM62',
                '--evalue', '0.001',
                '-o', args.alignment_file,
                '-p', str(args.nproc),
                '--block-size', str(args.block_size),
                '--index-chunks', str(args.index_chunks),
                '--tmpdir', args.tmpdir,
                '--compress', compression]

        if not args.verbose:
            command += ['--quiet']

        if args.sensitive:
            command += ['--sensitive']

        subprocess.check_call(command)
    except:
        message = 'DIAMOND finished abnormally.'
        give_user_feedback(message, args.log_file, args.quiet, error=True)

        sys.exit(1)

    if args.compress:
        setattr(args, 'alignment_file', '{0}.gz'.format(args.alignment_file))

    message = 'Homology search done! File {0} created.'.format(
            args.alignment_file)
    give_user_feedback(message, args.log_file, args.quiet)

    return


def import_contig_names(fasta_file, log_file, quiet):
    message = 'Importing contig names from {0}.'.format(fasta_file)
    give_user_feedback(message, log_file, quiet)
    
    contig_names = set()
    
    with open(fasta_file, 'r') as f1:
        for line in f1:
            if line.startswith('>'):
                contig = line.split(' ')[0].lstrip('>').rstrip()
                
                if contig in contig_names:
                    message = (
                            'it looks like your fasta file contains duplicate '
                            'headers! The first duplicate encountered is {0}, '
                            'but there might be more...'.format(contig))
                    give_user_feedback(message, log_file, quiet, error=True)
                    
                    sys.exit(1)
                    
                contig_names.add(contig)

    return contig_names


def import_ORFs(proteins_fasta, log_file, quiet):
    message = 'Parsing ORF file {0}'.format(proteins_fasta)
    give_user_feedback(message, log_file, quiet)

    contig2ORFs = {}
    
    with open(proteins_fasta, 'r') as f1:
        for line in f1:
            line = line.rstrip()

            if line.startswith('>'):
                ORF = line.split(' ')[0].lstrip('>')
                contig = ORF.rsplit('_', 1)[0]

                if contig not in contig2ORFs:
                    contig2ORFs[contig] = []

                contig2ORFs[contig].append(ORF)

    return contig2ORFs


def parse_tabular_alignment(
        alignment_file, one_minus_r, log_file, quiet):
    message = 'Parsing alignment file {0}.'.format(alignment_file)
    give_user_feedback(message, log_file, quiet)

    compressed = False
    if alignment_file.endswith('.gz'):
        compressed = True

        f1 = gzip.open(alignment_file, 'rb')
    else:
        f1 = open(alignment_file, 'r')

    ORF2hits = {}
    all_hits = set()

    ORF = 'first ORF'
    ORF_done = False
    for line in f1:
        if compressed:
            line = line.decode('utf-8')

        if line.startswith(ORF) and ORF_done == True:
            # The ORF has already surpassed its minimum allowed bit-score.
            continue

        line = line.rstrip().split('\t')

        if not line[0] == ORF:
            # A new ORF is reached.
            ORF = line[0]
            best_bitscore = decimal.Decimal(line[11])
            ORF2hits[ORF] = []

            ORF_done = False

        bitscore = decimal.Decimal(line[11])
        
        if bitscore >= one_minus_r * best_bitscore:
            # The hit has a high enough bit-score to be included.
            hit = line[1]

            ORF2hits[ORF].append((hit, bitscore),)
            all_hits.add(hit)
        else:
            # The hit is not included because its bit-score is too low.
            ORF_done = True

    f1.close()
                
    return (ORF2hits, all_hits)


if __name__ == '__main__':
    sys.exit('Run \'CAT\' to run CAT or BAT.')
