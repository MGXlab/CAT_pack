#!/usr/bin/env python3

import argparse
import os
import sys

import about
import check
import shared


def parse_arguments():
    parser = argparse.ArgumentParser(prog='CAT summarise',
                                     description='Summarise a named CAT or '
                                                 'BAT classification file.',
                                     usage='CAT summarise -i -o (-c) '
                                           '[options] [-h / --help]',
                                     add_help=False)
    
    required = parser.add_argument_group('Required arguments')
    
    required.add_argument('-i',
                          '--input_file',
                          dest='named_input_file',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to named CAT contig classification file '
                               'or BAT bin classification file. Currently '
                               'only official ranks are supported, and only '
                               'classification files containing a single '
                               'classification per contig / bin.')
    required.add_argument('-o',
                          '--output_file',
                          dest='output_file',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to output file.')
    
    optional = parser.add_argument_group('Optional arguments')
    
    optional.add_argument('-c',
                          '--contigs_fasta',
                          dest='contigs_fasta',
                          metavar='',
                          required=False,
                          type=str,
                          help='Path to contigs fasta file. This is required '
                               'if you want to summarise a contig '
                               'classification file.')
    optional.add_argument('--force',
                          dest='force',
                          required=False,
                          action='store_true',
                          help='Force overwrite existing files.')
    optional.add_argument('-q',
                          '--quiet',
                          dest='quiet',
                          required=False,
                          action='store_true',
                          help='Suppress verbosity.')
    optional.add_argument('-h',
                          '--help',
                          action='help',
                          help='Show this help message and exit.')

    (args, extra_args) = parser.parse_known_args()

    extra_args = [arg for (i, arg) in enumerate(extra_args) if
                  (i, arg) != (0, 'summarise')]
    if len(extra_args) > 0:
        sys.exit('error: too much arguments supplied:\n{0}'
                 ''.format('\n'.join(extra_args)))

    return args


def import_contig_lengths(contigs_fasta, log_file, quiet):
    message = 'Gathering contig lengths from {0}.'.format(contigs_fasta)
    shared.give_user_feedback(message, log_file, quiet)

    contig2length = {}

    with open(contigs_fasta, 'r') as f1:
        for line in f1:
            line = line.rstrip()

            if line.startswith('>'):
                contig = line.split(' ')[0].lstrip('>')

                contig2length[contig] = 0
            else:
                try:
                    contig2length[contig] += len(line)
                except:
                    message = ('ERROR: {0} is not a contigs fasta'
                               ''.format(contigs_fasta))
                    shared.give_user_feedback(message,
                                              log_file,
                                              quiet,
                                              error=True)

                    sys.exit(1)

    return contig2length


def summarise_contigs(input_file, output_file, contigs_fasta, force, quiet):
    # Currently summarise does not a allow for a log file.
    log_file = None
    
    message = '# CAT v{0}.'.format(about.__version__)
    shared.give_user_feedback(message, log_file, quiet, show_time=False)

    errors = []

    errors.append(check.check_input_file(input_file, log_file, quiet))

    if not force:
        errors.append(check.check_output_file(output_file, log_file, quiet))

    if True in errors:
        sys.exit(1)
        
    contig2length = import_contig_lengths(contigs_fasta, log_file, quiet)

    message = 'Summarising...'
    shared.give_user_feedback(message, log_file, quiet)

    with open(input_file, 'r') as f1:
        for line in f1:
            if line.startswith('#'):
                line = line.split('\t')
                
                if line[0] != '# contig':
                    message = ('ERROR: {0} is not a CAT classification file.'
                               ''.format(input_file))
                    shared.give_user_feedback(message,
                                              log_file,
                                              quiet,
                                              error=True)

                    if line[0] == '# bin':
                        message = ('ERROR: {0} appears to be a BAT '
                                   'classification file. If you want to '
                                   'summarise bin classifications, just '
                                   'don\'t supply a contigs fasta and '
                                   'everything should be fine!'
                                   ''.format(input_file))
                        shared.give_user_feedback(message,
                                                  log_file,
                                                  quiet,
                                                  error=True)
                        
                    sys.exit(1)
                    
                try:
                    superkingdom_index = line.index('superkingdom')
                except:
                    message = ('ERROR: official ranks not found in header of '
                               '{0}. Make sure that the CAT classification '
                               'file is named with official ranks with \'CAT '
                               'add_names --only_official\'.'
                               ''.format(input_file))
                    shared.give_user_feedback(message,
                                              log_file,
                                              quiet,
                                              error=True)

                    sys.exit(1)

                break
        else:
            message = 'ERROR: input file does not have a recognisable header.'
            shared.give_user_feedback(message, log_file, quiet, error=True)

            sys.exit(1)
            
    length = {}
    length['unclassified'] = []

    ORFs = {}

    official_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family',
                      'genus', 'species']

    for rank in official_ranks:
        length[rank] = {}
        ORFs[rank] = {}
    
    n = 0
    contig_trace = set()
    doubles = set()
    with open(input_file, 'r') as f1:
        for line in f1:
            line = line.rstrip()

            if line.startswith('#'):
                continue

            n += 1

            line = line.split('\t')

            contig = line[0]

            if contig in contig_trace:
                doubles.add(contig)

            contig_trace.add(contig)
            
            if contig not in contig2length:
                message = ('ERROR: contig {0} in CAT classification file is '
                           'not found in supplied contigs fasta file. Are you '
                           'sure the CAT classification file is based on the '
                           'contigs fasta?'.format(contig))
                shared.give_user_feedback(message, log_file, quiet, error=True)

                sys.exit(1)

            if line[1] == 'unclassified':
                length['unclassified'].append(contig2length[contig])

                continue

            for (i, classification) in enumerate(line[superkingdom_index:]):
                classification = classification.rsplit(': ', 1)[0].rstrip('*')
                
                rank = official_ranks[i]

                if classification not in length[rank]:
                    length[rank][classification] = []

                    ORFs[rank][classification] = []

                length[rank][classification].append(contig2length[contig])

                # Note that the total number of ORFs on a contig is reproted,
                # not only the number of ORFs a classification is based on.
                ORFs_on_contig = int(line[2].split('/')[1].split(' ')[0])
                ORFs[rank][classification].append(ORFs_on_contig)

    if len(doubles) != 0:
        message = ('ERROR: some contigs have multiple classifications. CAT '
                   'summarise currently does not allow for this. Contigs with '
                   'multiple classifications: {0}.'
                   ''.format(', '.join(list(doubles))))
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)
        
    if n != len(contig2length):
        message = ('ERROR: the number of classified contigs is not the same '
                   'as the number of contigs in contigs fasta. Are you sure '
                   'the CAT classification file is based on the contigs '
                   'fasta?')
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    with open(output_file, 'w') as outf1:
        number_of_contigs = len(contig2length)
        total_length = sum(contig2length.values())
        number_of_classified_contigs = number_of_contigs - len(length['unclassified'])
        total_classified_length = total_length - sum(length['unclassified'])

        outf1.write('# total number of contigs in {0} is {1} representing {2} '
                    'positions.\n'
                    ''.format(contigs_fasta,
                              number_of_contigs,
                              total_length))
        outf1.write('# {0} contigs are classified ({1:.2f}%) representing {2} '
                    'positions ({3:.2f}%) in {4}.\n'
                    ''.format(number_of_classified_contigs,
                              number_of_classified_contigs / number_of_contigs * 100,
                              total_classified_length,
                              total_classified_length / total_length * 100,
                              input_file))
        outf1.write('#\n')
        outf1.write('# rank\t'
                    'clade\t'
                    'number of contigs\t'
                    'number of ORFs\t'
                    'number of positions\n')

        for rank in official_ranks:
            for clade in sorted(length[rank],
                               key=lambda x: sum(length[rank][x]),
                               reverse=True):
                outf1.write('{0}\t{1}\t{2}\t{3}\t{4}\n'
                            ''.format(rank,
                                      clade,
                                      len(length[rank][clade]),
                                      sum(ORFs[rank][clade]),
                                      sum(length[rank][clade])))

    message = '{0} is created!'.format(output_file)
    shared.give_user_feedback(message, log_file, quiet)
    
    
def summarise_bins(input_file, output_file, force, quiet):
    # Currently summarise does not a allow for a log file.
    log_file = None
    
    message = '# CAT v{0}.'.format(about.__version__)
    shared.give_user_feedback(message, log_file, quiet, show_time=False)
    
    errors = []

    errors.append(check.check_input_file(input_file, log_file, quiet))

    if not force:
        errors.append(check.check_output_file(output_file, log_file, quiet))

    if True in errors:
        sys.exit(1)
        
    message = 'Summarising...'
    shared.give_user_feedback(message, log_file, quiet)

    with open(input_file, 'r') as f1:
        for line in f1:
            if line.startswith('#'):
                line = line.split('\t')
                
                if line[0] != '# bin':
                    message = ('ERROR: {0} is not a BAT classification file.'
                               ''.format(input_file))
                    shared.give_user_feedback(message,
                                              log_file,
                                              quiet,
                                              error=True)

                    if line[0] == '# contig':
                        message = ('ERROR: {0} appears to be a CAT '
                                   'classification file. If you want to '
                                   'summarise contig classifications, please '
                                   'supply a contigs fasta.'
                                   ''.format(input_file))
                        shared.give_user_feedback(message,
                                                  log_file,
                                                  quiet,
                                                  error=True)
                        
                    sys.exit(1)
                    
                try:
                    superkingdom_index = line.index('superkingdom')
                except:
                    message = ('ERROR: official ranks not found in header of '
                               '{0}. Make sure that the BAT classification '
                               'file is named with official ranks with \'CAT '
                               'add_names --only_official\'.'
                               ''.format(input_file))
                    shared.give_user_feedback(message,
                                              log_file,
                                              quiet,
                                              error=True)

                    sys.exit(1)

                break
        else:
            message = 'ERROR: input file does not have a recognisable header.'
            shared.give_user_feedback(message, log_file, quiet, error=True)

            sys.exit(1)
            
    number_of_bins = {}
    number_of_bins['unclassified'] = 0
    
    official_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family',
                      'genus', 'species']

    for rank in official_ranks:
        number_of_bins[rank] = {}
    
    n = 0
    bin_trace = set()
    doubles = set()
    with open(input_file, 'r') as f1:
        for line in f1:
            line = line.rstrip()

            if line.startswith('#'):
                continue

            n += 1

            line = line.split('\t')

            bin_ = line[0]

            if bin_ in bin_trace:
                doubles.add(bin_)

            bin_trace.add(bin_)
            
            if line[1] == 'unclassified':
                number_of_bins['unclassified'] += 1
                
                continue

            for (i, classification) in enumerate(line[superkingdom_index:]):
                classification = classification.rsplit(': ', 1)[0].rstrip('*')
                
                rank = official_ranks[i]

                if classification not in number_of_bins[rank]:
                    number_of_bins[rank][classification] = 0

                number_of_bins[rank][classification] += 1
                
    if len(doubles) != 0:
        message = ('ERROR: some bins have multiple classifications. CAT '
                   'summarise currently does not allow for this. Bins with '
                   'multiple classifications: {0}.'
                   ''.format(', '.join(list(doubles))))
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)
        
    number_of_classified_bins = n - number_of_bins['unclassified']

    with open(output_file, 'w') as outf1:
        outf1.write('# total number of bins is {0}, of which {1} ({2:.2f}%) '
                    'are classified.\n'
                    ''.format(n,
                              number_of_classified_bins,
                              number_of_classified_bins / n * 100))
        outf1.write('#\n')
        outf1.write('# rank\tclade\tnumber of bins\n')

        for rank in official_ranks:
            for clade in sorted(number_of_bins[rank],
                               key=lambda x: number_of_bins[rank][x],
                               reverse=True):
                outf1.write('{0}\t{1}\t{2}\n'
                            ''.format(rank,
                                      clade,
                                      number_of_bins[rank][clade]))
                            
    message = '{0} is created!'.format(output_file)
    shared.give_user_feedback(message, log_file, quiet)
    
    
def summarise(args):
    (input_file,
     output_file,
     contigs_fasta,
     force,
     quiet) = check.convert_arguments(args)
    
    if contigs_fasta == None:
        summarise_bins(input_file, output_file, force, quiet)
    else:
        summarise_contigs(input_file, output_file, contigs_fasta, force, quiet)
        
        
def run():
    args = parse_arguments()

    summarise(args)
    
    
if __name__ == '__main__':
    sys.exit('Please run \'CAT summarise\' to summarise a named CAT contig '
             'classification file or named BAT bin classification file.')
