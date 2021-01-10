#!/usr/bin/env/ python3

import hashlib
import os
import subprocess
import sys

import shared


def check_md5_gz(gz_file, md5_file, log_file, quiet):
    message = 'Checking file integrity via MD5 checksum.'
    shared.give_user_feedback(message, log_file, quiet)

    with open(md5_file, 'r') as f:
        md5_exp = f.read().split(' ')[0]

    if md5_exp == '':
        message = ('WARNING: no MD5 found in {0}. Integrity of {1} can not be '
                'established.'.format(md5_file, gz_file))
        shared.give_user_feedback(message, log_file, quiet)
    else:
        md5 = hashlib.md5()

        block_size = 4096
        with open(gz_file, 'rb') as f:
            for chunk in iter(lambda: f.read(block_size), b''):
                md5.update(chunk)
        md5 = md5.hexdigest()

        if md5 != md5_exp:
            message = 'MD5 of {0} does not check out.'.format(gz_file)
            shared.give_user_feedback(message, log_file, quiet, error=True)

            sys.exit(1)
        else:
            message = 'MD5 of {0} checks out.'.format(gz_file)
            shared.give_user_feedback(message, log_file, quiet)

    return


def check_memory(Gb):
    total_memory = None
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

    if os.path.isdir(out_prefix):
        message = 'prefix for output files ({0}) is a directory.'.format(
                out_prefix)
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    dir_ = out_prefix.rsplit('/', 1)[0]

    if not os.path.isdir(dir_):
        message = ('can not find output directory {0} to which output files '
                'should be written.'.format(dir_))
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
        message = ('can not find Prodigal. Please check whether it is '
                'installed or the path to the binaries is provided.')
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
        message = ('can not find DIAMOND. Please check whether it is '
                'installed or the path to the binaries is provided.')
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_bin_folder(bin_folder, bin_suffix, log_file, quiet):
    error = False

    if not os.path.isdir(bin_folder):
        message = 'can not find the bin folder.'
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

        return error
    
    tmp = []
    for file_ in os.listdir(bin_folder):
        if file_.startswith('.'):
            # Skip hidden files.
            continue

        if not file_.endswith(bin_suffix):
            continue

        if '.concatenated.' in file_:
            # Skip concatenated contig fasta and predicted protein fasta files
            # from earlier runs.
            continue
        
        tmp.append(file_)

    if len(tmp) == 0:
        message = (
                'no bins found with suffix {0} in bin folder. You can set the '
                'suffix with the [-s / --bin_suffix] argument.'.format(
                    bin_suffix))
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True
    elif len(tmp) == 1:
        message = (
                'WARNING: a single bin is found. You can run BAT in single '
                'bin mode, with \'CAT bin\' as opposed to \'CAT bins\' for a '
                'set of bins. Both modes will give the same results, but you '
                'might find single mode more convenient for your workflow.')
        shared.give_user_feedback(message, log_file, quiet)

    return error


def check_bin_fasta(bin_fasta, log_file, quiet):
    error = False

    if check_fasta(bin_fasta, log_file, quiet):
        error = True

    if os.path.isdir(bin_fasta):
        message = (
                '{0} is a directory. If you want to classify more than 1 bin '
                'you can run \'CAT bins\' instead of \'CAT bin\'.'.format(
                    bin_fasta))
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_folders_for_run(
        taxonomy_folder,
        nodes_dmp,
        names_dmp,
        database_folder,
        diamond_database,
        fastaid2LCAtaxid_file,
        taxids_with_multiple_offspring_file,
        step_list,
        log_file,
        quiet):
    error = False

    if not os.path.isdir(taxonomy_folder):
        message = 'can not find the taxonomy folder.'
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True
    else:
        if not nodes_dmp or not names_dmp:
            message = ('nodes.dmp and / or names.dmp not found in the '
                    'taxonomy folder.')
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

    if not os.path.isdir(database_folder):
        message = 'can not find the database folder.'
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True
    else:
        if not diamond_database and 'align' in step_list:
            message = 'DIAMOND database not found in database folder.'
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

        if not fastaid2LCAtaxid_file:
            message = 'file fastaid2LCAtaxid is not found in database folder.'
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

        if not taxids_with_multiple_offspring_file:
            message = ('file taxids_with_multiple_offspring not found in '
                    'database folder.')
            shared.give_user_feedback(message, log_file, quiet, error=True)

            error = True

    return error


def check_output_file(output_file, log_file, quiet):
    error = False

    if os.path.isfile(output_file):
        message = (
                'output file {0} already exists. You can choose to overwrite '
                'existing files with the [--force] argument.'.format(
                    output_file))
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_input_file(input_file, log_file, quiet):
    error = False

    if not os.path.isfile(input_file):
        message = 'input file {0} does not exist.'.format(input_file)
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_in_and_output_file(input_file, output_file, log_file, quiet):
    error = False

    if input_file == output_file:
        message = 'input file and output file can not be the same.'
        shared.give_user_feedback(message, log_file, quiet, error=True)

        error = True

    return error


def check_top(top, r, log_file, quiet):
    error = False

    if top < 50:
        message = (
                'WARNING: [--top] is set lower than 50. This might conflict '
                'with future runs with higher settings of the '
                '[-r / --range] parameter, see README.md.')
        shared.give_user_feedback(message, log_file, quiet)
        
    if top <= r:
        message = '[--top] should be higher than [-r / --range].'
        shared.give_user_feedback(message, log_file, quiet, error=True)
        
        error = True

    return error


def check_fasta(file_, log_file, quiet):
    error = False

    if not os.path.isfile(file_):
        error = True
    else:
        with open(file_, 'r') as f1:
            for n, line in enumerate(f1):
                if n == 0:
                    if not line.startswith('>'):
                        error = True

                    break

    if error:
        message = '{0} is not a fasta file.'.format(file_)
        shared.give_user_feedback(message, log_file, quiet, error=True)

    return error


def check_whether_ORFs_are_based_on_contigs(
        contig_names, contig2ORFs, log_file, quiet):
    for contig in contig2ORFs:
        if contig not in contig_names:
            message = (
                    'found a protein in the predicted proteins fasta file '
                    'that can not be traced back to one of the contigs in the '
                    'contigs fasta file: {0}. Proteins should be named '
                    'contig_name_#.'.format(contig2ORFs[contig][0]))
            shared.give_user_feedback(message, log_file, quiet, error=True)

            sys.exit(1)

    return
            
            
if __name__ == '__main__':
    sys.exit('Run \'CAT\' to run CAT or BAT.')
