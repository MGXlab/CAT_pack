#!/usr/bin/env python3

import argparse
import os
import sys

import about
import check
import shared


def parse_arguments():
    parser = argparse.ArgumentParser(prog='CAT summarise',
                                     description='Summarise a named CAT '
                                                 'contig classification file.',
                                     usage='CAT summarise -c -i -o '
                                           '[-h / --help]',
                                     add_help=False)
    
    required = parser.add_argument_group('Required arguments')
    
    required.add_argument('-c',
                          '--contigs_fasta',
                          dest='contigs_fasta',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to contigs fasta file.')
    required.add_argument('-i',
                          '--input_file',
                          dest='named_input_file',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to named CAT contig classification file. '
                               'Currently only official ranks are supported, '
                               'and only classification files containing a '
                               'single classification per contig.')
    required.add_argument('-o',
                          '--output_file',
                          dest='output_file',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to output file.')
    
    optional = parser.add_argument_group('Optional arguments')
    
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

    with open(contigs_fasta, 'r') as f:
        for line in f:
            line = line.rstrip()

            if line.startswith('>'):
                contig = line.split(' ')[0].lstrip('>')

                contig2length[contig] = 0
            else:
                contig2length[contig] += len(line)

    return contig2length


def summarise(args):
    (contigs_fasta,
     input_file,
     output_file,
     quiet) = check.convert_arguments(args)
    
    # Currently summarise does not a allow for a log file.
    log_file = None
    
    message = '# CAT v{0}.'.format(about.__version__)
    shared.give_user_feedback(message, log_file, quiet, show_time=False)
    
    if not os.path.isfile(input_file):
        message = 'ERROR: input file does not exist.'
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    contig2length = import_contig_lengths(contigs_fasta, log_file, quiet)

    message = 'Summarising...'
    shared.give_user_feedback(message, log_file, quiet)

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                line = line.split('\t')

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
    with open(input_file, 'r') as f:
        for line in f:
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

            if line[1].startswith('unclassified'):
                length['unclassified'].append(contig2length[contig])

                continue

            for (i, classification) in enumerate(line[superkingdom_index:]):
                classification = classification.rsplit(': ', 1)[0].rstrip('*')
                
                rank = official_ranks[i]

                if classification not in length[rank]:
                    length[rank][classification] = []

                    ORFs[rank][classification] = []

                length[rank][classification].append(contig2length[contig])

                ORFs[rank][classification].append(int(line[2]))

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

    with open(output_file, 'w') as outf:
        number_of_contigs = len(contig2length)
        total_length = sum(contig2length.values())
        number_of_classified_contigs = number_of_contigs - len(length['unclassified'])
        total_classified_length = total_length - sum(length['unclassified'])

        outf.write('# total number of contigs in {0} is {1} representing {2} '
                   'positions.\n'
                   ''.format(contigs_fasta,
                             number_of_contigs,
                             total_length))
        outf.write('# {0} contigs are classified ({1:.2f}%) representing {2} '
                   'positions ({3:.2f}%) in {4}.\n'
                   ''.format(number_of_classified_contigs,
                             number_of_classified_contigs / number_of_contigs * 100,
                             total_classified_length,
                             total_classified_length / total_length * 100,
                             input_file))
        outf.write('#\n')
        outf.write('# classification level\t'
                   'name\t'
                   'number of contigs\t'
                   'number of ORFs\t'
                   'number of positions\n')

        for rank in official_ranks:
            for name in sorted(length[rank],
                               key=lambda x: sum(length[rank][x]),
                               reverse=True):
                outf.write('{0}\t{1}\t{2}\t{3}\t{4}\n'
                           ''.format(rank,
                                     name,
                                     len(length[rank][name]),
                                     sum(ORFs[rank][name]),
                                     sum(length[rank][name])))

    message = '{0} is created!'.format(output_file)
    shared.give_user_feedback(message, log_file, quiet)


def run():
    args = parse_arguments()

    summarise(args)


if __name__ == '__main__':
    sys.exit('Please run \'CAT summarise\' to summarise a named CAT contig '
             'classification file.')
