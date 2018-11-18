#!/usr/bin/env python3

import argparse
import datetime
import multiprocessing
import os
import sys

import about
import check
import shared
import tax


def parse_arguments():
    parser = argparse.ArgumentParser(prog='CAT bins',
                                     description='Run Bin Annotation Tool '
                                                 '(BAT).',
                                     usage='CAT bins -b -d -t '
                                           '[options] [-h / --help]',
                                     add_help=False)
    
    required = parser.add_argument_group('Required arguments')

    required.add_argument('-b',
                          '--bin_folder',
                          dest='bin_folder',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to directory containing bins.')
    required.add_argument('-d',
                          '--database_folder',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to folder that contains database files.')
    required.add_argument('-t',
                          '--taxonomy_folder',
                          dest='taxonomy_folder',
                          metavar='',
                          required=True,
                          type=str,
                          help='Path to folder that contains taxonomy files.')

    optional = parser.add_argument_group('Optional arguments')
    
    optional.add_argument('-s',
                          '--bin_suffix',
                          dest='bin_suffix',
                          metavar='',
                          required=False,
                          type=str,
                          default='.fna',
                          help='Suffix of bins in bin folder (default: .fna).')
    optional.add_argument('-r',
                          '--range',
                          dest='r',
                          metavar='',
                          required=False,
                          type=float,
                          choices = [i for i in range(51)],
                          default=5,
                          help='r parameter [0-50] (default: 5).')
    optional.add_argument('-f',
                          '--fraction',
                          dest='f',
                          metavar='',
                          required=False,
                          type=float,
                          choices = [i / 100 for i in range(0, 100)],
                          default=0.3,
                          help='f parameter [0-0.99] (default: 0.3).')
    optional.add_argument('-o',
                          '--out_prefix',
                          dest='out_prefix',
                          metavar='',
                          required=False,
                          type=str,
                          default='out.BAT',
                          help='Prefix for output files (default: out.BAT).')
    optional.add_argument('-p',
                          '--proteins_fasta',
                          dest='predicted_proteins_fasta',
                          metavar='',
                          required=False,
                          type=str,
                          help='Path to concatenated predicted proteins fasta '
                               'file generated during an earlier run of BAT. '
                               'If supplied, BAT will skip the protein '
                               'prediction step.')
    optional.add_argument('-a',
                          '--diamond_alignment',
                          dest='diamond_file',
                          metavar='',
                          required=False,
                          type=str,
                          help='Path to Diamond alignment table generated '
                               'during an earlier run of BAT. If supplied, '
                               'BAT will skip the Diamond alignment step and '
                               'only classify the bins. A concatenated '
                               'predicted proteins fasta file should also be '
                               'supplied with argument [-p / --proteins].')
    optional.add_argument('-n',
                          '--nproc',
                          dest='nproc',
                          metavar='',
                          required=False,
                          type=int,
                          default=multiprocessing.cpu_count(),
                          help='Number of cores to deploy by Diamond '
                               '(default: maximum).')
    optional.add_argument('--path_to_prodigal',
                          dest='path_to_prodigal',
                          metavar='',
                          required=False,
                          type=str,
                          default='prodigal',
                          help='Path to Prodigal binaries. Please supply if '
                               'BAT can not find Prodigal.')
    optional.add_argument('--path_to_diamond',
                          dest='path_to_diamond',
                          metavar='',
                          required=False,
                          type=str,
                          default='diamond',
                          help='Path to Diamond binaries. Please supply if '
                               'BAT can not find Diamond.')
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

    (args, extra_args) = parser.parse_known_args()

    extra_args = [arg for (i, arg) in enumerate(extra_args) if
                  (i, arg) != (0, 'bins')]
    if len(extra_args) > 0:
        sys.exit('error: too much argument supplied:\n{0}'
                 ''.format('\n'.join(extra_args)))

    return args


def import_bins(bin_folder,
                bin_suffix,
                log_file,
                quiet):
    message = 'Importing bins from {0}/.'.format(bin_folder)
    shared.give_user_feedback(message, log_file, quiet)

    bin2contigs = {}
    contig_names = set()

    for file in os.listdir(bin_folder):
        if not file.endswith(bin_suffix):
            continue

        if '.concatenated.fa' in file:
            # This serves to skip the concatenated contig fasta and predicted
            # protein fasta from earlier runs if BAT is run with -p and -a, and
            # the bin folder is the same folder as the folder in which CAT its
            # executables are located, and -s is .faa or .fasta.
            continue

        bin_ = file.rsplit(bin_suffix, 1)[0]

        bin2contigs[bin_] = []

        with open('{0}/{1}'.format(bin_folder, file), 'r') as f1:
            for line in f1:
                if line.startswith('>'):
                    contig = line.split(' ')[0].rstrip().lstrip('>')

                    # Add bin name in front of the contig name.
                    new_contig_name = '{0}_{1}'.format(bin_, contig)
                    
                    if new_contig_name in contig_names:
                        message = ('ERROR: BAT has encountered {0} twice in '
                                   'bin {1}. Each fasta header should be '
                                   'unique in each bin.'
                                   ''.format(contig, bin_))
                        shared.give_user_feedback(message,
                                                  log_file,
                                                  quiet,
                                                  error=True)

                        sys.exit(1)

                    contig_names.add(new_contig_name)

                    bin2contigs[bin_].append(new_contig_name)
                    
    message = '{0} bins found!'.format(len(bin2contigs))
    shared.give_user_feedback(message, log_file, quiet)

    return (bin2contigs, contig_names)


def make_concatenated_fasta(concatenated_fasta,
                            bin_folder,
                            bin_suffix,
                            log_file,
                            quiet):
    message = 'Writing {0}.'.format(concatenated_fasta)
    shared.give_user_feedback(message, log_file, quiet)

    files = os.listdir(bin_folder)

    with open(concatenated_fasta, 'w') as outf1:
        for file in files:
            if not file.endswith(bin_suffix):
                continue

            bin_ = file.rsplit(bin_suffix, 1)[0]

            with open('{0}/{1}'.format(bin_folder, file), 'r') as f1:
                for line in f1:
                    if line.startswith('>'):
                        contig = line.split(' ')[0].rstrip().lstrip('>')
                        
                        # Add bin name in front of the contig name.
                        outf1.write('>{0}_{1}\n'.format(bin_, contig))
                    else:
                        outf1.write(line)
                        
                        
def bins(args):
    step_list = []

    (bin_folder,
     database_folder,
     taxonomy_folder,
     bin_suffix,
     one_minus_r,
     f,
     out_prefix,
     predicted_proteins_fasta,
     diamond_file,
     nproc,
     path_to_prodigal,
     path_to_diamond,
     quiet,
     no_log) = check.convert_arguments(args)
    
    if no_log:
        log_file = None
    else:
        log_file = '{0}.log'.format(out_prefix)
        with open(log_file, 'w') as outf1:
            pass
        
    message = '# CAT v{0}.'.format(about.__version__)
    shared.give_user_feedback(message, log_file, quiet, show_time=False)
    
    # Check at which state to start.
    if predicted_proteins_fasta is None and diamond_file is None:
        message = ('\n'
                   'BAT is running. Protein prediction, alignment, and bin '
                   'classification are carried out.\n'
                   'Rarw!\n\n'
                   'Supplied command: {0}\n\n'
                   'Bin folder: {1}/\n'
                   'Taxonomy folder: {2}/\n'
                   'Database folder: {3}/\n'
                   'Parameter r: {4}\n'
                   'Parameter f: {5}\n'
                   'Log file: {6}\n\n'
                   '-----------------\n'.format(' '.join(sys.argv),
                                                bin_folder,
                                                taxonomy_folder,
                                                database_folder,
                                                args.r,
                                                args.f,
                                                log_file))
        shared.give_user_feedback(message, log_file, quiet, show_time=False)

        step_list.append('run_prodigal')
        step_list.append('run_diamond')
    elif (predicted_proteins_fasta is not None
          and diamond_file is None):
        message = ('\n'
                   'BAT is running. Since a predicted protein fasta is '
                   'supplied, only alignment and bin classification are '
                   'carried out.\n'
                   'Rarw!\n\n'
                   'Supplied command: {0}\n\n'
                   'Bin folder: {1}/\n'
                   'Taxonomy folder: {2}/\n'
                   'Database folder: {3}/\n'
                   'Parameter r: {4}\n'
                   'Parameter f: {5}\n'
                   'Log file: {5}\n\n'
                   '-----------------\n'.format(' '.join(sys.argv),
                                                bin_folder,
                                                taxonomy_folder,
                                                database_folder,
                                                args.r,
                                                args.f,
                                                log_file))
        shared.give_user_feedback(message, log_file, quiet, show_time=False)

        step_list.append('run_diamond')
    elif (predicted_proteins_fasta is not None and
          diamond_file is not None):
        message = ('\n'
                   'BAT is running. Since a predicted protein fasta and '
                   'Diamond alignment file are supplied, only bin '
                   'classification is carried out.\n'
                   'Rarw!\n\n'
                   'Supplied command: {0}\n\n'
                   'Bin folder: {1}/\n'
                   'Taxonomy folder: {2}/\n'
                   'Database folder: {3}/\n'
                   'Parameter r: {4}\n'
                   'Parameter f: {5}\n'
                   'Log file: {6}\n\n'
                   '-----------------\n'.format(' '.join(sys.argv),
                                                bin_folder,
                                                taxonomy_folder,
                                                database_folder,
                                                args.r,
                                                args.f,
                                                log_file))
        shared.give_user_feedback(message, log_file, quiet, show_time=False)
    elif (predicted_proteins_fasta is None and
          diamond_file is not None):
        message = ('ERROR: if you want BAT to only classify a set of bins, '
                   'you should not only supply a Diamond alignment table but '
                   'also a concatenated predicted protein fasta file with '
                   'argument [-p / --proteins].')
        shared.give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)
        
    # Do binaries, taxonomy folder, and database folder check and 
    # set parameters.
    message = 'Doing some pre-flight checks first.'
    shared.give_user_feedback(message, log_file, quiet, show_time=False)

    errors = []

    errors.append(check.check_bin_folder(bin_folder,
                                         bin_suffix,
                                         log_file,
                                         quiet))
    
    errors.append(check.check_out_prefix(out_prefix, log_file, quiet))
    
    if 'run_prodigal' in step_list:
        errors.append(check.check_prodigal_binaries(path_to_prodigal,
                                                    log_file,
                                                    quiet))
        concatenated_fasta = '{0}.concatenated.fasta'.format(out_prefix)
        predicted_proteins_fasta = ('{0}.concatenated.predicted_proteins.faa'
                                    ''.format(out_prefix))
        predicted_proteins_gff = ('{0}.concatenated.predicted_proteins.gff'
                                  ''.format(out_prefix))
        
    if 'run_diamond' in step_list:
        errors.append(check.check_diamond_binaries(path_to_diamond,
                                                   log_file,
                                                   quiet))
        diamond_file = ('{0}.concatenated.alignment.diamond'
                        ''.format(out_prefix))
    else:
        diamond_file = diamond_file
        
    errors.append(check.check_folders_for_run(taxonomy_folder,
                                              database_folder,
                                              step_list,
                                              log_file,
                                              quiet))
    
    bin2classification_output_file = ('{0}.bin2classification.txt'
                                      ''.format(out_prefix))
    ORF2LCA_output_file = '{0}.ORF2LCA.txt'.format(out_prefix)

    errors.append(check.check_output_files(bin2classification_output_file,
                                           ORF2LCA_output_file,
                                           log_file,
                                           quiet))

    if 'run_prodigal' not in step_list:
        if not check.check_if_file_is_fasta(predicted_proteins_fasta):
            message = ('ERROR: {0} is not a fasta file.'
                       ''.format(predicted_proteins_fasta))
            shared.give_user_feedback(message, log_file, quiet, error=True)

            errors.append(True)
            
    if True in errors:
        sys.exit(1)
        
    (nodes_dmp,
     names_dmp,
     prot_accession2taxid_file) = check.inspect_taxonomy_folder(taxonomy_folder)
    (nr_file,
     diamond_database,
     fastaid2LCAtaxid_file,
     taxids_with_multiple_offspring_file) = check.inspect_database_folder(database_folder)

    message = 'Ready to fly!\n\n-----------------\n'
    shared.give_user_feedback(message, log_file, quiet, show_time=False)
    
    # Start BAT.
    (bin2contigs, contig_names) = import_bins(bin_folder,
                                              bin_suffix,
                                              log_file,
                                              quiet)

    if 'run_prodigal' in step_list:
        make_concatenated_fasta(concatenated_fasta,
                                bin_folder,
                                bin_suffix,
                                log_file,
                                quiet)

        shared.run_prodigal(path_to_prodigal,
                            concatenated_fasta,
                            predicted_proteins_fasta,
                            predicted_proteins_gff,
                            log_file,
                            quiet)

    if 'run_diamond' in step_list:
        shared.run_diamond(path_to_diamond,
                           diamond_database,
                           predicted_proteins_fasta,
                           diamond_file,
                           nproc,
                           log_file,
                           quiet)
        
    contig2ORFs = shared.import_ORFs(predicted_proteins_fasta,
                                     log_file,
                                     quiet)
    
    check.check_whether_ORFs_are_based_on_contigs(contig_names,
                                                  contig2ORFs,
                                                  log_file,
                                                  quiet)
    
    (ORF2hits,
     all_hits) = shared.parse_diamond_file(diamond_file,
                                           one_minus_r,
                                           log_file,
                                           quiet)

    (taxid2parent, taxid2rank) = tax.import_nodes(nodes_dmp, log_file, quiet)
    fastaid2LCAtaxid = tax.import_fastaid2LCAtaxid(fastaid2LCAtaxid_file,
                                                   all_hits,
                                                   log_file,
                                                   quiet)
    taxids_with_multiple_offspring = tax.import_taxids_with_multiple_offspring(taxids_with_multiple_offspring_file,
                                                                               log_file,
                                                                               quiet)
    
    message = ('BAT is flying! Files {0} and {1} are created.'
               ''.format(bin2classification_output_file,
                         ORF2LCA_output_file))
    shared.give_user_feedback(message, log_file, quiet)

    with open(bin2classification_output_file, 'w') as outf1, open(ORF2LCA_output_file, 'w') as outf2:
        outf1.write('# bin\tclassification\tnumber of ORFs in bin\t'
                    'number of ORFs classification is based on\tlineage\t'
                    'lineage scores\n')
        outf2.write('# ORF\tbin\tlineage\tbit-score\n')
        
        for bin_ in sorted(bin2contigs):
            LCAs_ORFs = []

            for contig in sorted(bin2contigs[bin_]):
                if contig not in contig2ORFs:
                    continue

                for ORF in contig2ORFs[contig]:
                    if ORF not in ORF2hits:
                        outf2.write('{0}\t{1}\tORF has no hit to database.\n'
                                    ''.format(ORF, bin_))

                        continue

                    (taxid,
                     top_bitscore) = tax.find_LCA_for_ORF(ORF2hits[ORF],
                                                          fastaid2LCAtaxid,
                                                          taxid2parent)
                     
                    if taxid.startswith('no taxid found'):
                        outf2.write('{0}\t{1}\t{2}\t{3}\n'
                                    ''.format(ORF,
                                              bin_,
                                              taxid,
                                              top_bitscore))
                    else:
                        lineage = tax.find_lineage(taxid, taxid2parent)
                        starred_lineage = tax.star_lineage(lineage,
                                                           taxids_with_multiple_offspring)

                        outf2.write('{0}\t{1}\t{2}\t{3}\n'
                                    ''.format(ORF,
                                              bin_,
                                              ';'.join(starred_lineage[::-1]),
                                              top_bitscore))
                                       
                    LCAs_ORFs.append((taxid, top_bitscore))
                    
            if len(LCAs_ORFs) == 0:
                outf1.write('{0}\tunclassified (no hits to database)\n'
                            ''.format(bin_))

                continue

            (lineages,
             lineages_scores,
             based_on_number_of_ORFs) = tax.find_weighted_LCA(LCAs_ORFs,
                                                              taxid2parent,
                                                              f)

            if lineages == 'no ORFs with taxids found.':
                outf1.write('{0}\tunclassified '
                            '(hits not found in taxonomy files)\n'
                            ''.format(bin_))

                continue

            if lineages == 'no lineage whitelisted.':
                outf1.write('{0}\tunclassified '
                            '(no lineage reached minimum bit-score support)\n'
                            ''.format(bin_))

                continue
            
            # The bin has a valid classification.
            total_number_of_ORFs = sum([len(contig2ORFs[contig]) for
                                        contig in bin2contigs[bin_] if
                                        contig in contig2ORFs])
            
            for (i, lineage) in enumerate(lineages):
                starred_lineage = tax.star_lineage(lineage,
                                                   taxids_with_multiple_offspring)
                
                scores = ['{0:.2f}'.format(score) for score in
                          lineages_scores[i]]
                
                if len(lineages) == 1:
                    # There is only one classification.
                    outf1.write('{0}\tclassified\t{1}\t{2}\t{3}\t{4}\n'
                                ''.format(bin_,
                                          total_number_of_ORFs,
                                          based_on_number_of_ORFs,
                                          ';'.join(starred_lineage[::-1]),
                                          ';'.join(scores[::-1])))
                else:
                    # There are multiple classifications.
                    outf1.write('{0}\tclassified ({1}/{2})'
                                '\t{3}\t{4}\t{5}\t{6}\n'
                                ''.format(bin_,
                                          i + 1, len(lineages),
                                          total_number_of_ORFs,
                                          based_on_number_of_ORFs,
                                          ';'.join(starred_lineage[::-1]),
                                          ';'.join(scores[::-1])))
                                   
    message = ('\n-----------------\n'
               '[{0}] BAT is done! {1} bins classified.'
               ''.format(datetime.datetime.now(), len(bin2contigs)))
    shared.give_user_feedback(message, log_file, quiet, show_time=False)
  
    if f < 0.5:
        message = ('WARNING: since f is set to smaller than 0.5, one bin '
                   'may have multiple classifications.')
        shared.give_user_feedback(message, log_file, quiet, show_time=False)
        
        
def run():
    args = parse_arguments()

    bins(args)


if __name__ == '__main__':
    sys.exit('Please run \'CAT bins\' to run Bin Annotation Tool (BAT).')
