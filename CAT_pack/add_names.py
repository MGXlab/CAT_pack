#!/usr/bin/env python3

import argparse
import os
import sys

import about
import check
import shared
import tax


def parse_arguments():
    parser = argparse.ArgumentParser(prog='CAT add_names',
                                     description='Add taxonomic names to CAT '
                                                 'or BAT output files.',
                                     usage='CAT add_names -i -o -t '
                                           '[options] [-h / --help]',
                                     add_help=False)
    
    required = parser.add_argument_group('Required arguments')
    
    required.add_argument('-i',
                          '--input_file',
                          dest='input_file',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to input file. Can be either '
                               'classification output file or ORF2LCA output '
                               'file.')
    required.add_argument('-o',
                          '--output_file',
                          dest='output_file',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to output file.')
    required.add_argument('-t',
                          '--taxonomy_folder',
                          dest='taxonomy_folder',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to folder that contains taxonomy files.')
    
    optional = parser.add_argument_group('Optional arguments')
    
    optional.add_argument('--only_official',
                          dest='only_official',
                          required=False,
                          action='store_true',
                          help='Only output official level names.')
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
                  (i, arg) != (0, 'add_names')]
    if len(extra_args) > 0:
        sys.exit('error: too much arguments supplied:\n{0}'
                 ''.format('\n'.join(extra_args)))

    return args


def add_names(args):
    (input_file,
     output_file,
     taxonomy_folder,
     only_official,
     force,
     quiet) = check.convert_arguments(args)

    # Currently add_names does not allow for a log file.
    log_file = None
    
    message = '# CAT v{0}.'.format(about.__version__)
    shared.give_user_feedback(message, log_file, quiet, show_time=False)

    errors = []

    errors.append(check.check_input_file(input_file, log_file, quiet))

    if not force:
        errors.append(check.check_output_file(output_file, log_file, quiet))

    if True in errors:
        sys.exit(1)
        
    (nodes_dmp,
     names_dmp,
     prot_accession2taxid_file) = check.inspect_taxonomy_folder(taxonomy_folder)

    (taxid2parent, taxid2rank) = tax.import_nodes(nodes_dmp, log_file, quiet)
    taxid2name = tax.import_names(names_dmp, log_file, quiet)

    message = 'Appending names...'
    shared.give_user_feedback(message, log_file, quiet)

    with open(input_file, 'r') as f1:
        for line in f1:
            if line.startswith('#'):
                line = line.rstrip().split('\t')

                try:
                    lineage_index = line.index('lineage')
                except:
                    message = ('ERROR: {0} is not a supported classification '
                               'file.'.format(input_file))
                    shared.give_user_feedback(message,
                                              log_file,
                                              quiet,
                                              error=True)

                    sys.exit(1)
                    
                try:
                    scores_index = line.index('lineage scores')
                except:
                    scores_index = None

                full_length = len(line)

                break
        else:
            message = ('ERROR: {0} is not a supported classification file.'
                       ''.format(input_file))
            shared.give_user_feedback(message, log_file, quiet, error=True)

            sys.exit(1)
            
    with open(input_file, 'r') as f1, open(output_file, 'w') as outf1:
        for line in f1:
            line = line.rstrip()

            if line.startswith('#'):
                if only_official:
                    outf1.write('{0}\tsuperkingdom\tphylum\tclass\torder\t'
                                'family\tgenus\tspecies\n'.format(line))
                else:
                    outf1.write('{0}\tfull lineage names\n'.format(line))
                    
                continue
            
            line = line.split('\t')

            if len(line) != full_length:
                # Entry does not have a full annotation.
                outf1.write('{0}\n'.format('\t'.join(line)))

                continue

            if (line[1].startswith('no taxid found') or
                line[2].startswith('no taxid found')):
                # ORF has database hits but the accession number is not found
                # in the taxonomy files.
                outf1.write('{0}\n'.format('\t'.join(line)))

                continue
            
            lineage = line[lineage_index].split(';')

            if scores_index:
                scores = line[scores_index].split(';')
            else:
                scores = None

            if only_official:
                names = tax.convert_to_official_names(lineage,
                                                      taxid2rank,
                                                      taxid2name,
                                                      scores)
            else:
                names = tax.convert_to_names(lineage,
                                             taxid2rank,
                                             taxid2name,
                                             scores)

            outf1.write('{0}\t{1}\n'.format('\t'.join(line), '\t'.join(names)))

    message = 'Names written to {0}!'.format(output_file)
    shared.give_user_feedback(message, log_file, quiet)
    
    
def run():
    args = parse_arguments()

    add_names(args)
    
    
if __name__ == '__main__':
    sys.exit('Please run \'CAT add_names\' to add taxonomic names to CAT or '
             'BAT output files.')
