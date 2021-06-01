#!/usr/bin/env python3

import argparse
import sys

import about
import check
import shared


def parse_arguments():
    parser = argparse.ArgumentParser(
            prog='CAT summarise',
            description='Summarise a named CAT or BAT classification file.',
            usage='CAT summarise -i -o (-c) [options] [-h / --help]',
            add_help=False)
    
    required = parser.add_argument_group('Required arguments')
    shared.add_argument(required, 'input_file', True,
            help_=(
                'Path to named CAT contig classification file or BAT bin '
                'classification file. Currently only official ranks are '
                'supported, and only classification files containing a single '
                'classification per contig / bin. If you want to summarise a '
                'contig classification file, you have to supply the contigs '
                'fasta file with argument [-c / --contigs_fasta].'))
    shared.add_argument(required, 'output_file', True)
    
    optional = parser.add_argument_group('Optional arguments')
    shared.add_argument(optional, 'contigs_fasta', False,
            help_=('Path to contigs fasta file. Required if you want to '
                'summarise a contig classification file.'))
    shared.add_argument(optional, 'force', False)
    shared.add_argument(optional, 'quiet', False)
    shared.add_argument(optional, 'help', False)

    (args, extra_args) = parser.parse_known_args()

    extra_args = [arg for (i, arg) in enumerate(extra_args) if
                  (i, arg) != (0, 'summarise')]
    if len(extra_args) > 0:
        sys.exit('error: too much arguments supplied:\n{0}'.format(
            '\n'.join(extra_args)))

    # Add extra arguments.
    shared.expand_arguments(args)

    return args


def import_contig_lengths(contigs_fasta, log_file, quiet):
    message = 'Gathering contig lengths from {0}.'.format(contigs_fasta)
    shared.give_user_feedback(message, log_file, quiet)

    contig2length = {}

    with open(contigs_fasta, 'r') as f1:
        for line in f1:
            line = line.rstrip()

            if line.startswith('>'):
                contig = line.split()[0].lstrip('>')

                contig2length[contig] = 0
            else:
                try:
                    contig2length[contig] += len(line)
                except:
                    message = '{0} is not a contigs fasta'.format(
                            contigs_fasta)
                    shared.give_user_feedback(
                            message, log_file, quiet, error=True)

                    sys.exit(1)

    return contig2length


def summarise_contigs(args):
    message = '# CAT v{0}.'.format(about.__version__)
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)

    errors = []

    errors.append(
            check.check_input_file(args.input_file, args.log_file, args.quiet))

    if not args.force:
        errors.append(
                check.check_output_file(
                    args.output_file, args.log_file, args.quiet))

    errors.append(
            check.check_in_and_output_file(
                args.input_file, args.output_file, args.log_file, args.quiet))

    if True in errors:
        sys.exit(1)
        
    contig2length = import_contig_lengths(
            args.contigs_fasta, args.log_file, args.quiet)

    message = 'Summarising...'
    shared.give_user_feedback(message, args.log_file, args.quiet)

    with open(args.input_file, 'r') as f1:
        for line in f1:
            if line.startswith('#'):
                line = line.split('\t')
                
                if line[0] != '# contig':
                    message = '{0} is not a CAT classification file.'.format(
                            args.input_file)
                    shared.give_user_feedback(
                            message, args.log_file, args.quiet, error=True)

                    if line[0] == '# bin':
                        message = (
                                '{0} appears to be a BAT classification file. '
                                'If you want to summarise bin '
                                'classifications, simply don\'t supply a '
                                'contigs fasta and everything should be fine.'
                                ''.format(args.input_file))
                        shared.give_user_feedback(
                                message, args.log_file, args.quiet, error=True)
                        
                    sys.exit(1)
                    
                try:
                    superkingdom_index = line.index('superkingdom')
                except:
                    message = (
                            'official ranks not found in header of {0}. Make '
                            'sure that the CAT classification file is named '
                            'with official ranks with \'CAT add_names '
                            '--only_official\'.'.format(args.input_file))
                    shared.give_user_feedback(
                            message, args.log_file, args.quiet, error=True)

                    sys.exit(1)

                break
        else:
            message = 'input file does not have a recognisable header.'
            shared.give_user_feedback(message, args.log_file, args.quiet,
                    error=True)

            sys.exit(1)
            
    length = {}
    length['no taxid assigned'] = []

    ORFs = {}

    official_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family',
                      'genus', 'species']

    for rank in official_ranks:
        length[rank] = {}
        ORFs[rank] = {}
    
    n = 0
    contig_trace = set()
    doubles = set()
    with open(args.input_file, 'r') as f1:
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
                message = (
                        'contig {0} in CAT classification file is not found '
                        'in supplied contigs fasta file. Are you sure the CAT '
                        'classification file is based on the contigs fasta?'
                        ''.format(contig))
                shared.give_user_feedback(message, args.log_file, args.quiet,
                        error=True)

                sys.exit(1)

            if line[1] == 'no taxid assigned':
                length['no taxid assigned'].append(contig2length[contig])

                continue

            for (i, classification) in enumerate(line[superkingdom_index:]):
                classification = classification.rsplit(': ', 1)[0].rstrip('*')
                
                rank = official_ranks[i]

                if classification not in length[rank]:
                    length[rank][classification] = []

                    ORFs[rank][classification] = []

                length[rank][classification].append(contig2length[contig])

                # NOTE that the total number of ORFs on a contig is reported,
                # not only the number of ORFs a classification is based on.
                ORFs_on_contig = int(line[2].split('/')[1].split(' ')[0])
                ORFs[rank][classification].append(ORFs_on_contig)

    if len(doubles) != 0:
        message = (
                'some contigs have multiple classifications. CAT summarise '
                'currently does not allow for this. Contigs with multiple '
                'classifications: {0}.'.format(', '.join(list(doubles))))
        shared.give_user_feedback(message, args.log_file, args.quiet,
                error=True)

        sys.exit(1)
        
    if n != len(contig2length):
        message = (
                'the number of classified contigs is not the same as the '
                'number of contigs in contigs fasta. Are you sure the CAT '
                'classification file is based on the contigs fasta?')
        shared.give_user_feedback(message, args.log_file, args.quiet,
                error=True)

        sys.exit(1)

    with open(args.output_file, 'w') as outf1:
        n_contigs = len(contig2length)
        total_length = sum(contig2length.values())
        n_classified_contigs = n_contigs - len(length['no taxid assigned'])
        total_classified_length = total_length - sum(
                length['no taxid assigned'])

        outf1.write('# total number of contigs in {0} is {1:,d} representing '
                '{2:,d} positions.\n'.format(
                    args.contigs_fasta, n_contigs, total_length))
        outf1.write('# {0:,d} contigs have taxonomy assigned ({1:.2f}%) '
                'representing {2:,d} positions ({3:.2f}%) in {4}.\n'.format(
                    n_classified_contigs,
                    n_classified_contigs / n_contigs * 100,
                    total_classified_length,
                    total_classified_length / total_length * 100,
                    args.input_file))
        outf1.write('#\n')
        outf1.write(
                '# rank\t'
                'clade\t'
                'number of contigs\t'
                'number of ORFs\t'
                'number of positions\n')

        for rank in official_ranks:
            for clade in sorted(length[rank],
                    key=lambda x: sum(length[rank][x]), reverse=True):
                outf1.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(
                    rank,
                    clade,
                    len(length[rank][clade]),
                    sum(ORFs[rank][clade]),
                    sum(length[rank][clade])))

    message = '{0} is created!'.format(args.output_file)
    shared.give_user_feedback(message, args.log_file, args.quiet)

    return
    
    
def summarise_bins(args):
    message = '# CAT v{0}.'.format(about.__version__)
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)
    
    errors = []

    errors.append(
            check.check_input_file(args.input_file, args.log_file, args.quiet))

    if not args.force:
        errors.append(
                check.check_output_file(
                    args.output_file, args.log_file, args.quiet))

    errors.append(
            check.check_in_and_output_file(
                args.input_file, args.output_file, args.log_file, args.quiet))

    if True in errors:
        sys.exit(1)
        
    message = 'Summarising...'
    shared.give_user_feedback(message, args.log_file, args.quiet)

    with open(args.input_file, 'r') as f1:
        for line in f1:
            if line.startswith('#'):
                line = line.split('\t')
                
                if line[0] != '# bin':
                    message = '{0} is not a BAT classification file.'.format(
                            args.input_file)
                    shared.give_user_feedback(
                            message, args.log_file, args.quiet, error=True)

                    if line[0] == '# contig':
                        message = (
                                '{0} appears to be a CAT classification file. '
                                'If you want to summarise contig '
                                'classifications, supply a contigs fasta with '
                                'argument [-c / --contigs_fasta].'.format(
                                    args.input_file))
                        shared.give_user_feedback(
                                message, args.log_file, args.quiet, error=True)
                        
                    sys.exit(1)
                    
                try:
                    superkingdom_index = line.index('superkingdom')
                except:
                    message = (
                            'official ranks not found in header of {0}. Make '
                            'sure that the BAT classification file is named '
                            'with official ranks with \'CAT add_names '
                            '--only_official\'.'.format(args.input_file))
                    shared.give_user_feedback(
                            message, args.log_file, args.quiet, error=True)

                    sys.exit(1)

                break
        else:
            message = 'input file does not have a recognisable header.'
            shared.give_user_feedback(message, args.log_file, args.quiet,
                    error=True)

            sys.exit(1)
            
    n_bins = {}
    n_bins['no taxid assigned'] = 0
    
    official_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family',
                      'genus', 'species']

    for rank in official_ranks:
        n_bins[rank] = {}
    
    n = 0
    bin_trace = set()
    doubles = set()
    with open(args.input_file, 'r') as f1:
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
            
            if line[1] == 'no taxid assigned':
                n_bins['no taxid assigned'] += 1
                
                continue

            for (i, classification) in enumerate(line[superkingdom_index:]):
                classification = classification.rsplit(': ', 1)[0].rstrip('*')
                
                rank = official_ranks[i]

                if classification not in n_bins[rank]:
                    n_bins[rank][classification] = 0

                n_bins[rank][classification] += 1
                
    if len(doubles) != 0:
        message = (
                'some bins have multiple classifications. CAT summarise '
                'currently does not allow for this. Bins with multiple '
                'classifications: {0}.'.format(', '.join(list(doubles))))
        shared.give_user_feedback(message, args.log_file, args.quiet,
                error=True)

        sys.exit(1)
        
    n_classified_bins = n - n_bins['no taxid assigned']

    with open(args.output_file, 'w') as outf1:
        outf1.write('# total number of bins is {0:,d}, of which {1:,d} '
                '({2:.2f}%) have taxonomy assigned.\n'.format(
                    n, n_classified_bins, n_classified_bins / n * 100))
        outf1.write('#\n')
        outf1.write('# rank\tclade\tnumber of bins\n')

        for rank in official_ranks:
            for clade in sorted(n_bins[rank],
                    key=lambda x: n_bins[rank][x], reverse=True):
                outf1.write('{0}\t{1}\t{2}\n'.format(
                    rank, clade, n_bins[rank][clade]))
                            
    message = '{0} is created!'.format(args.output_file)
    shared.give_user_feedback(message, args.log_file, args.quiet)

    return
    
    
def run():
    args = parse_arguments()

    if not args.contigs_fasta:
        summarise_bins(args)
    else:
        summarise_contigs(args)

    return


if __name__ == '__main__':
    sys.exit('Run \'CAT summarise\' to summarise a named CAT contig '
            'classification file or named BAT bin classification file.')
