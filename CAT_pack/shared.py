#!/usr/bin/env python3

import datetime
import math
import subprocess
import sys


def give_user_feedback(message,
                       log_file=None,
                       quiet=False,
                       show_time=True,
                       error=False):
    time = datetime.datetime.now()

    if show_time:
        message = '[{0}] {1}\n'.format(time, message)
    else:
        message = '{0}\n'.format(message)

    if log_file:
        with open(log_file, 'a') as outf:
            outf.write(message)

    if not quiet and not error:
        sys.stdout.write(message)

    if not quiet and error:
        sys.stderr.write(message)
        
        
def run_prodigal(path_to_prodigal,
                 contigs_fasta,
                 predicted_proteins_fasta,
                 predicted_proteins_gff,
                 log_file,
                 quiet):
    message = ('Running Prodigal for ORF prediction. Files {0} and {1} will '
               'be generated. Do not forget to cite Prodigal when using CAT '
               'or BAT in your publication!'.format(predicted_proteins_fasta,
                                                    predicted_proteins_gff))
    give_user_feedback(message, log_file, quiet)

    try:
        subprocess.check_call([path_to_prodigal,
                               '-i', contigs_fasta,
                               '-a', predicted_proteins_fasta,
                               '-o', predicted_proteins_gff,
                               '-p', 'meta',
                               '-q',
                               '-f', 'gff'])
    except:
        message = 'ERROR: Prodigal finished abnormally.'
        give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    message = 'ORF prediction done!'
    give_user_feedback(message, log_file, quiet)
    
    
def run_diamond(path_to_diamond,
                diamond_database,
                predicted_proteins_fasta,
                diamond_file,
                nproc,
                log_file,
                quiet):
    message = ('Homology search with Diamond is starting. Please be patient. '
               'Do not forget to cite Diamond when using CAT or BAT in your '
               'publication!\n'
               '\t\t\t\tquery: {0}\n'
               '\t\t\t\tdatabase: {1}'.format(predicted_proteins_fasta,
                                              diamond_database))
    give_user_feedback(message, log_file, quiet)

    try:
        subprocess.check_call([path_to_diamond,
                               'blastp',
                               '-d', diamond_database,
                               '-q', predicted_proteins_fasta,
                               '--top', '50',
                               '-o', diamond_file,
                               '-p', str(nproc),
                               '--quiet'])
    except:
        message = 'ERROR: Diamond finished abnormally.'
        give_user_feedback(message, log_file, quiet, error=True)

        sys.exit(1)

    message = 'Homology search done! File {0} created.'.format(diamond_file)
    give_user_feedback(message, log_file, quiet)


def import_ORFs(predicted_proteins_fasta, log_file, quiet):
    message = 'Parsing ORF file {0}'.format(predicted_proteins_fasta)
    give_user_feedback(message, log_file, quiet)

    contig2ORFs = {}
    
    with open(predicted_proteins_fasta, 'r') as f:
        for line in f:
            line = line.rstrip()

            if line.startswith('>'):
                ORF = line.split(' ')[0].lstrip('>')
                contig = ORF.rsplit('_', 1)[0]

                if contig not in contig2ORFs:
                    contig2ORFs[contig] = []

                contig2ORFs[contig].append(ORF)

    return contig2ORFs


def parse_diamond_file(diamond_file,
                       bitscore_fraction_cutoff_1,
                       log_file,
                       quiet):
    message = 'Parsing Diamond file {0}.'.format(diamond_file)
    give_user_feedback(message, log_file, quiet)

    ORF2hits = {}
    all_hits = set()

    ORF = 'first ORF'
    ORF_done = False
    with open(diamond_file, 'r') as f:
        for line in f:
            if line.startswith(ORF) and ORF_done == True:
                # The ORF has already surpassed its minimum allowed bitscore.
                continue

            line = line.rstrip().split('\t')

            if not line[0] == ORF:
                # A new ORF is reached.
                ORF = line[0]
                best_bitscore = float(line[11])
                ORF2hits[ORF] = []

                ORF_done = False

            bitscore = float(line[11])
            
            if bitscore >= bitscore_fraction_cutoff_1 * best_bitscore:
                # The hit has a high enough bitscore to be included.
                hit = line[1]

                ORF2hits[ORF].append((hit, bitscore))
                all_hits.add(hit)
            else:
                # The hit is not included because its bitscore is too low.
                ORF_done = True
                
    return (ORF2hits, all_hits)


if __name__ == '__main__':
    sys.exit('Please run \'CAT\' to run CAT or BAT.')
