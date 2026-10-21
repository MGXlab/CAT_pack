#!/usr/bin/env python3

import subprocess
import sys
from dataclasses import field, dataclass
from pathlib import Path

from CAT_pack7 import tax
from CAT_pack7.config.settings import RatSettings
from CAT_pack7.tools.tool_runner import tool_output
from CAT_pack7.utils.errors import ExternalToolError
from CAT_pack7.utils.logging import Report


def process_CAT_table(settings: RatSettings, output_table: Path)-> dict[str, list[tuple[str,...]]]:
    """
    takes CAT contig2classification or bin2classification table and taxonomy
    nodes file as input; returns a dictionary that contains the complete
    lineage for each contig/bin:
    (list of tuples of the format [(tax_id1,tax_id2...),
                                   (rank1,rank2,...),
                                   (score1,score2...)])
    bin id contains the file extension!
    """
    c2c: dict[str, list[tuple[str, ...]]] = {}
    # get rank for each taxid
    taxid2rank=tax.read_nodes(settings.files.nodes, ranks=True)

    with open(output_table) as contigs:

        # @Tina, I moved the first line check out of the for loop and if the
        # first line has a # (which it has) it skips it with next()
        if not next(contigs).startswith("#"):
            contigs.seek(0) # go back to line 0 if there is no # found

        for contig in contigs:

            contig_cols = contig.strip().split('\t')

            # get contig id, lineage and lineage scores
            contig_id=contig_cols[0]
            c2c[contig_id]=[]

            # Build a list of tuples for the lineage:
            # (tax_id1,tax_id2...), (rank1,rank2,...), (score1,score2...)]
            if len(contig_cols)>3:
                tax_ids=contig_cols[3].strip().split(';')
                scores=contig_cols[4].strip().split(';')
                c2c[contig_id]=[tuple([t for t in tax_ids]),
                               tuple([taxid2rank[t.strip('*')] for t in tax_ids]),
                               tuple([s for s in scores])]

    return c2c



@dataclass
class MappingResult:
    read_dict: dict[str, dict[str, list[str]]] = field(default_factory=dict)
    unmapped_reads: dict[str, set[str]] = field(default_factory=lambda: {'fw': set(), 'rev': set()})
    sum_of_reads: int = 0
    paired: bool =False
    contig_dict: dict[str, int] = field(default_factory=dict)
    max_primary: int = 1

def process_bam_file(settings: RatSettings, report: Report) -> MappingResult:
    """
    # function that processes the bam file(s) by picking out the primary alignments
    # and storing the contigs that the forward and reverse reads map to. If only
    # one file is given, checks whether it's paired or not before it counts the
    # unmapped reads.
    # Keeps track of unmapped reads: If paired read file, then if checks the flag
    # for whether the unmapped read is the first or the second in the pair; if
    # it's single reads, it just fills up the fw key of the dictionary, and if
    # there is a reverse file, then the reverse file only adds to the rev key.
    """
    map_result = MappingResult()
    for index, bam_file in enumerate(settings.files.bam_files):
        # open alignment file with samtools, so that we don't need pysam dependency
        cmd=[str(settings.samtools), 'view', '-F', '0x900', bam_file]
        with tool_output("samtools") as output, subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=output, text=True) as proc:

            for mapping in proc.stdout:
                mapping=mapping.rstrip().split('\t')
                read_id, flag, contig, score=mapping[0], int(mapping[1]), mapping[2], int(mapping[4])
        
                if flag & 0x1:
                    map_result.paired=True

                map_result.read_dict.setdefault(read_id, {
                    "contig": [],
                    "bin": [],
                    "taxon_bin": [],
                    "taxon_contig": [],
                    "taxon_contig_dm": [],
                    "taxon_read_dm": [],
                })

                if flag & 0x4 or score < settings.mapping_quality:
                    direction = "rev" if index == 1 or flag & 0x80 else "fw"
                    map_result.unmapped_reads[direction].add(read_id)
                else:
                    map_result.read_dict[read_id]["contig"].append(contig)

            if proc.wait() != 0:
                raise ExternalToolError("samtools", f"view exited with a non 0 status: {proc.returncode}")

    contig2bin = invert_bin_dict(settings.files.bin2contigs)

    make_contig_dict(map_result, settings, contig2bin)


    if map_result.paired or len(settings.files.bam_files)>1:
        map_result.sum_of_reads=2*len(map_result.read_dict)
    else:
        map_result.sum_of_reads = len(map_result.read_dict)

    return map_result

# @Bastiaan and @Tina are inplace operations too complex? Return values can also be used.
def make_contig_dict(map_result: MappingResult, settings: RatSettings, contig2bin: dict[str, str]) -> None:
    max_primary = 2 if len(settings.files.bam_files) > 1 or map_result.paired else 1
    contig_dict = {"*": 0}

    for read in map_result.read_dict.values():
        contigs = read["contig"]
        contig_dict["*"] += max_primary - len(contigs)

        for contig in contigs:
            contig_dict[contig] = contig_dict.get(contig, 0) + 1

            if contig2bin:
                read["bin"].append(contig2bin.get(contig, "unbinned"))

    map_result.contig_dict = contig_dict
    map_result.max_primary = max_primary


#
# def make_contig_dict(read_dict,
#                      bam_files,
#                      paired,
#                      contig2bin):
#     """
#     takes the read dictionary that stores which read maps where and returns
#     a dictionary that stores how many reads bind to each contig. It also adds
#     bin annotation to the read dict
#     """
#     contig_dict={}
#     contig_dict['*']=0
#
#     # check how long the list of contigs per read is allowed to be to calculate
#     # unmapped reads and check that everything is working correctly
#     if len(bam_files)>1 or paired:
#         max_primary=2
#     else:
#         max_primary=1
#
#     for read in read_dict:
#         # if a read (pair) has multiple primary alignments, exit
#         if len(read_dict[read]['contig'])>max_primary:
#             message = (
#                     'Something is wrong with the alignment. Too many mappings '
#                     'detected. Please check your alignment file or map your reads '
#                     'using RAT.')
#             shared.give_user_feedback(message)
#             sys.exit(1)
#
#         # otherwise, determine the number of unmapped reads and add it to unmapped
#         # build a dictionary of contigs and how many reads map to them
#         # Add bin info to read dictionary
#         else:
#             contig_dict['*']+=max_primary-len(read_dict[read]['contig'])
#             for contig in read_dict[read]['contig']:
#                 if not contig in contig_dict:
#                     contig_dict[contig]=0
#
#                 contig_dict[contig]+=1
#                 if contig2bin:
#                     if contig in contig2bin:
#                         read_dict[read]['bin'].append(contig2bin[contig])
#                     else:
#                         read_dict[read]['bin'].append('unbinned')
#
#     return contig_dict, read_dict, max_primary
#
#
#
#
#
# def make_unclassified_seq_fasta(seq_fasta, unclassified_seq_ids,
#                                 unclassified_seq_fasta, f_format, f_mode,
#                                 suffix=''):
#
#
#     fasta_dict={}
#     if f_format == 'fasta':
#         with shared.optionally_compressed_handle(seq_fasta, 'r') as f1:
#             for line in f1:
#                 if line.startswith('>'):
#                     header = line.rstrip().split(' ')[0].lstrip('>')
#                     if header in fasta_dict:
#                         sys.exit('Duplicate fasta headers in the file!') # This should be prettified.
#                     fasta_dict.setdefault(header, '')
#                 else:
#                     fasta_dict[header] += line.rstrip()
#
#
#     elif f_format=='fastq':
#         with shared.optionally_compressed_handle(seq_fasta, 'r') as f1:
#             for n, line in enumerate(f1):
#                 if n % 4 == 0:
#                     if not line.startswith('@'):
#                         sys.exit('Unknown seuqence identifier symbol in line {n}!') # This should be prettified.
#                     header = line.rstrip().split(' ')[0].lstrip('@')
#
#                     if header in fasta_dict:
#                         sys.exit('Duplicate fasta headers in the file!')
#
#                     fasta_dict.setdefault(header, '')
#                 elif n % 4 == 1:
#                     fasta_dict[header] += line.rstrip()
#     else:
#         sys.exit('Unknown read file format!')
#
#     with open(unclassified_seq_fasta, f_mode) as outf:
#         # print(unclassified_seq_ids)
#         suffices=['/1', '/2', '_1', '_2']
#         for seq in unclassified_seq_ids:
#             if seq in fasta_dict:
#                 seq_id=seq
#             else:
#                 for s in suffices:
#                     if seq+s in fasta_dict:
#                         seq_id=seq+s
#         # print(seq)
#             outf.write('>{0}{1}\n{2}\n'.format(seq_id, suffix,
#                                                fasta_dict[seq]))
#
#
#
# def get_unclassified_contigs(contig2bin, c2c, b2c):
#     """
#     If the contig is in a bin that is unclassified, or if the contig is not in
#     a bin AND unclassified: put it in a set.
#     """
#     unclassified_contigs=list()
#     for c in c2c:
#         classified=False
#         if c in contig2bin:
#             if b2c[contig2bin[c]]:
#                 mag_lineage=b2c[contig2bin[c]][0]
#                 if len(mag_lineage)>2 or (len(mag_lineage)==2 and
#                                           mag_lineage[1]!='131567'):
#                     classified=True
#         if not classified:
#             if c2c[c]:
#                 lineage=c2c[c][0]
#                 if len(lineage)>2 or (len(lineage)==2 and
#                                           lineage[1]!='131567'):
#                     classified=True
#         if not classified:
#             unclassified_contigs.append(c)
#     return unclassified_contigs
#
#
#
# def get_contig_lengths(contig_file):
#     """
#     gets the lengths of all contigs from the contig file. Works for all types
#     of contig IDs, as long as it is a proper file with unique headers, because
#     it checks the length of the strings.
#     """
#     contig_length_dict={}
#     sum_of_nucleotides=0
#
#     # split the fasta file into contigs and prevent empty strings
#     # @Tina: opening a file like this will load the entire file into memory. I would open it line by line myself. This may get super slow depending on how you implement it, it can also be very fast.
#     with open(contig_file, 'r') as inp:
#         # @Tina: what does the [1:] do here? Also shouldn't the strip() remove the trailing newline? or do you mean rstrip()?
#         contigs=inp.read()[1:].strip().split('\n>')
#
#     # for each contig, store contig name and the length of sequence without whitespace
#     # for each contig, add its length to the sum of nucleotides in the assembly
#     for c in contigs:
#         c_id, c_seq=c.split('\n')[0].strip(), c.split('\n', 1)[1].strip().replace('\n','')
#         contig_length_dict[c_id]=len(c_seq)
#         sum_of_nucleotides+=len(c_seq)
#
#     return contig_length_dict, sum_of_nucleotides
#
#
#
#
#
# def process_bin_folder(path_to_folder, bin_suffix):
#     """
#     takes the path to a folder containing metagenomic bins and picks out all
#     nucleotide fasta files. Returns a dictionary that contains the bin ids as
#     keys and a set of all contig ids that belong in that bin as values.
#     Important: Only contains the binned sequences. Every contig that is not
#     in the output dictionary is unbinned.
#     @Tina: Don't we have a function from the BAT workflow that does this?
#     """
#     bin_dict={}
#     bins=os.listdir(path_to_folder)
#
#     for b in bins:
#         # grabs all fasta files and gets the bin id
#         if b.split('.')[-1]==bin_suffix.split('.')[-1] and not b.rsplit('.', 1)[0].endswith('unbinned'):
#             bin_id=b.strip()
#             bin_dict[bin_id]=set()
#
#             # open fasta file and split it into contigs
#             contigs=open(path_to_folder+b).read().split('>')[1:]
#
#             # add all contig ids to set
#
#             for contig in contigs:
#                 try:
#                     bin_dict[bin_id].add(contig.split()[0])
#                 except IndexError:
#                     sys.exit(f'The bin {b} is empty or not formatted correctly'
#                              ' (Look for a ">" sign at the end of the file). '
#                              'Did something go wrong during binning?')
#
#     return bin_dict


def invert_bin_dict(bin_dict):
    # WORKS

    """
    @Tina: this should work as well:
    d_rev = {v: k for k, v in d.items()}
    return d

    It's shorter, but of course you should make an extra check to see if the
    values of the first dictionary are unique. I do like unversal functions (e.g.
    invert_dict as opposed to invert_bin_dict.
    """
    
    contig2bin={}
    for b in bin_dict:
        for contig in bin_dict[b]:
            contig2bin[contig]=b
    return contig2bin


if __name__ == '__main__':
    sys.exit('Run \'CAT_pack reads\' to run Read Annotation Tool (RAT).')
