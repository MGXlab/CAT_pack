#!/usr/bin/env python3

import subprocess
import sys
from dataclasses import field, dataclass
from pathlib import Path
from typing import Iterator, Literal

from CAT_pack import tax
from CAT_pack.config.settings import RatSettings, ClassificationSettings
from CAT_pack.io.parsers import FastaParser, FastqParser
from CAT_pack.tools.aligner import run_aligner
from CAT_pack.tools.tool_runner import tool_output
from CAT_pack.utils.errors import ExternalToolError
from CAT_pack.utils.logging import Report
from CAT_pack.classification import ClassificationEngine
from CAT_pack.io.parsers import optionally_compressed_handle, AlignmentParser
from CAT_pack.io.writers import ClassificationWriter
from CAT_pack.utils.typing import EntityType


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

    with optionally_compressed_handle(output_table, 'r') as contigs:

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

    def iter_unmapped_read_ids(self) -> Iterator[str]:
        for read in self.unmapped_reads["fw"]:
            yield f"{read}_1"

        for read in self.unmapped_reads["rev"]:
            yield f"{read}_2"

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


def make_unclassified_seq_fasta(settings: RatSettings, map_result: MappingResult,
                                unclassified_contigs: list[str]) -> None:
    unclassified_seq_ids = set(unclassified_contigs)

    with settings.files.unclassified_fasta.open("w") as outf:
        for record in FastaParser(settings.files.contigs_fasta):
            if record.name in unclassified_seq_ids:
                outf.write(f">{record.name}\n{record.sequence}\n")

        for index, seq_fasta in enumerate(settings.files.read_files, start=1):
            direction = "fw" if index == 1 else "rev"
            unmapped_reads = map_result.unmapped_reads[direction]

            for record in FastqParser(seq_fasta):
                seq_id = record.name

                if seq_id not in unmapped_reads and seq_id.endswith(("/1", "/2", "_1", "_2")):
                    seq_id = seq_id[:-2]

                if seq_id in unmapped_reads:
                    outf.write(f">{seq_id}_{index}\n{record.sequence}\n")


def get_unclassified_contigs(contig2bin, c2c, b2c):
    """
    If the contig is in a bin that is unclassified, or if the contig is not in
    a bin AND unclassified: put it in a set.
    """
    unclassified_contigs=list()
    for c in c2c:
        classified=False
        if c in contig2bin:
            if b2c[contig2bin[c]]:
                mag_lineage=b2c[contig2bin[c]][0]
                if len(mag_lineage)>2 or (len(mag_lineage)==2 and
                                          mag_lineage[1]!='131567'):
                    classified=True
        if not classified:
            if c2c[c]:
                lineage=c2c[c][0]
                if len(lineage)>2 or (len(lineage)==2 and
                                          lineage[1]!='131567'):
                    classified=True
        if not classified:
            unclassified_contigs.append(c)
    return unclassified_contigs


def get_classification(*classifications):
    for classification in classifications:
        if classification:
            lineage = classification[0]
            if len(lineage) > 2 or (
                len(lineage) == 2
                and lineage[1].rstrip("*") != "131567"
            ):
                return classification
    return []

def get_contig_lengths(contig_file: Path) -> tuple[dict[str, int], int]:
    """
    gets the lengths of all contigs from the contig file. Works for all types
    of contig IDs, as long as it is a proper file with unique headers, because
    it checks the length of the strings.
    """
    contig_lengths = {contig.name: len(contig.sequence) for contig in FastaParser(contig_file)}
    return contig_lengths, sum(contig_lengths.values())
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

@dataclass
class Abundance:
    n_reads: int = 0
    n_nucleotides: int = 0
    classification: list[tuple[str, ...]] = field(default_factory=list)


@dataclass
class RatResult:
    mapping: MappingResult
    c2c: dict[str, list[tuple[str, ...]]]
    b2c: dict[str, list[tuple[str, ...]]]
    unclassified_contigs: list[str]
    all_unclassified: list[str]
    contig2bin: dict[str, str] = field(default_factory=dict)
    u2c: dict[str, list[tuple[str, ...]]] = field(default_factory=dict)
    contig_data: dict[str, Abundance] = field(default_factory=dict)
    taxon_dict: dict[str, Abundance] = field(default_factory=dict)
    bin_data: dict[str, Abundance] = field(default_factory=dict)


def load_classifications(settings: RatSettings, report: Report, map_result: MappingResult) -> RatResult:
    files = settings.files
    b2c = {}

    if files.contigs2classification is not None:
        c2c = process_CAT_table(settings, files.contigs2classification)
    else:
        c2c = {name: [] for name in FastaParser(files.contigs_fasta).iter_headers()}
    report.progress_next()

    if files.bin2classification is not None:
        b2c = process_CAT_table(settings, files.bin2classification)
    report.progress_next()

    contig2bin = invert_bin_dict(files.bin2contigs)
    unclassified_contigs = get_unclassified_contigs(
        contig2bin, c2c, b2c,
    )
    report.progress_next()

    all_unclassified = unclassified_contigs.copy()
    all_unclassified.extend(map_result.iter_unmapped_read_ids())
    report.progress_next()

    return RatResult(map_result, c2c, b2c, unclassified_contigs, all_unclassified, contig2bin)


def extract_unclassified_sequences(settings: RatSettings, report: Report, result: RatResult) -> None:
    make_unclassified_seq_fasta(
        settings, result.mapping, result.unclassified_contigs,
    )


def run_rat_alignment(settings: RatSettings, report: Report) -> None:
    run_aligner(settings.aligner, report)

# TODO: maybe integrate this with the run_classification function in pipeline.py
def write_unmapped2classification(settings: RatSettings, report: Report, result: RatResult) -> None:
    files = settings.files
    alignment = AlignmentParser(files.alignment_unmapped, settings.range_, workers=settings.threads).parse()
    report.progress_next()

    taxonomy = tax.Taxonomy.import_nodes(files.nodes)
    fastaid2LCAtaxid = tax.import_fastaid2LCAtaxid(files.fastaid2LCA, set(alignment.hit_names), workers=settings.threads)
    taxids_with_multiple_offspring = tax.import_taxids_with_multiple_offspring(files.import_taxids_with_multiple_offspring)
    report.progress_next()

    classifications = ClassificationEngine(
        taxonomy,
        taxonomy.taxid_of_hits(alignment.hit_names, fastaid2LCAtaxid),
        settings.fraction, settings.threads,
    ).classify(
        EntityType.SEQUENCE,
        {unclassified: [unclassified] for unclassified in result.all_unclassified},
        alignment
    )
    report.progress_next()

    with files.unmapped2classification.open("w") as unmapped_out:
        writer = ClassificationWriter(
            unmapped_out,
            entity_type=EntityType.SEQUENCE,
            fraction=settings.fraction,
            range_=settings.range_,
            branches=taxids_with_multiple_offspring,
            no_stars=settings.no_stars
        )
        for entity in range(len(classifications)):
            writer.write(classifications, entity)

    report.progress_next()
    return

def calculate_abundances(
    settings: RatSettings, report: Report, rat_result: RatResult,
) -> None:
    files = settings.files
    if files.unmapped2classification is not None:
        rat_result.u2c = process_CAT_table(
            settings, files.unmapped2classification,
        )
    report.progress_next()

    contig_length_dict, _ = get_contig_lengths(files.contigs_fasta)
    rat_result.contig_data = {
        contig: Abundance(n_reads, contig_length_dict[contig])
        for contig, n_reads in rat_result.mapping.contig_dict.items()
        if contig != "*"
    }
    rat_result.taxon_dict = {
        "unmapped": Abundance(rat_result.mapping.contig_dict["*"]),
        "unclassified": Abundance(),
    }
    rat_result.bin_data = {
        "unmapped": Abundance(rat_result.mapping.contig_dict["*"]),
    }
    report.progress_next()

    for contig, contig_data in rat_result.contig_data.items():
        bin_id = rat_result.contig2bin.get(contig, "unbinned")
        contig_data.classification = get_classification(
            rat_result.b2c.get(bin_id, []),
            rat_result.c2c.get(contig, []),
            rat_result.u2c.get(contig, []),
        )
        taxon = (
            ";".join(contig_data.classification[0])
            if contig_data.classification else "unclassified"
        )
        taxon_data = rat_result.taxon_dict.setdefault(
            taxon,
            Abundance(classification=contig_data.classification),
        )
        taxon_data.n_reads += contig_data.n_reads
        taxon_data.n_nucleotides += contig_data.n_nucleotides

        bin_data = rat_result.bin_data.setdefault(bin_id, Abundance())
        bin_data.n_reads += contig_data.n_reads
        bin_data.n_nucleotides += contig_data.n_nucleotides

    report.progress_next()

    for read in rat_result.mapping.iter_unmapped_read_ids():
        classification = get_classification(rat_result.u2c.get(read, []))
        if classification:
            taxon = ";".join(classification[0])
            taxon_data = rat_result.taxon_dict.setdefault(
                taxon, Abundance(classification=classification),
            )
            taxon_data.n_reads += 1
            rat_result.taxon_dict["unmapped"].n_reads -= 1

    report.progress_next()


def classify_reads(settings: RatSettings, report: Report, rat_result: RatResult) -> None:
    for name, read in rat_result.mapping.read_dict.items():
        read["taxon_bin"] = []
        read["taxon_contig"] = []
        read["taxon_contig_dm"] = []
        read["taxon_read_dm"] = []

        for contig, bin_id in zip(read["contig"], read["bin"]):
            for key, table, identifier in (
                ("taxon_bin", rat_result.b2c, bin_id),
                ("taxon_contig", rat_result.c2c, contig),
                ("taxon_contig_dm", rat_result.u2c, contig),
            ):
                classification = get_classification(
                    table.get(identifier, []),
                )
                read[key].append(
                    ";".join(classification[0]) if classification else "",
                )

        for number, direction in ((1, "fw"), (2, "rev")):
            if name in rat_result.mapping.unmapped_reads[direction]:
                classification = get_classification(
                    rat_result.u2c.get(f"{name}_{number}", []),
                )
                read["taxon_read_dm"].append(
                    f"{direction}: " + ";".join(classification[0])
                    if classification else "",
                )


def get_proper_fraction(rank_dict: dict[str, Abundance]) -> float:
    return sum(
        read.n_reads / read.n_nucleotides
        for read in rank_dict.values()
        if read.n_nucleotides
    )


def make_tax_table(settings: RatSettings, rat_result: RatResult) -> None:
    total_reads = rat_result.mapping.sum_of_reads
    divisor = get_proper_fraction(rat_result.contig_data)

    with settings.files.contig_abundance.open(
        "w", encoding="utf-8",
    ) as outf_contig:
        outf_contig.write(f"## command: {' '.join(sys.argv)}\n")
        outf_contig.write('# contig\tnumber of reads\tfraction of reads\t'
                          'contig_length\tAverage Coverage\tlineage\t'
                          'lineage ranks\tlineage scores\n')
        unmapped = rat_result.mapping.contig_dict["*"]
        read_frac = unmapped / total_reads if total_reads else 0.0
        outf_contig.write(
            f"unmapped\t{unmapped}\t{read_frac}\t0\t0\t\t\t\n",
        )

        for contig, read in rat_result.contig_data.items():
            read_frac = read.n_reads / total_reads if total_reads else 0.0
            corrected_frac = (
                read.n_reads / read.n_nucleotides / divisor
                if read.n_nucleotides and divisor else 0.0
            )
            lineage, ranks, scores = (
                (";".join(values) for values in read.classification)
                if read.classification else ("", "", "")
            )
            outf_contig.write(
                f"{contig}\t{read.n_reads}\t{read_frac}\t"
                f"{read.n_nucleotides}\t{corrected_frac}\t"
                f"{lineage}\t{ranks}\t{scores}\n"
            )

    with settings.files.complete_abundance.open("w", encoding="utf-8") as outf_complete:
        outf_complete.write(f"## command: {' '.join(sys.argv)}\n")
        outf_complete.write('# lineage\tnumber of reads\tfraction of reads\t'
                            'taxon length\tlineage ranks\n')

        for lineage, read in rat_result.taxon_dict.items():
            read_frac = read.n_reads / total_reads if total_reads else 0.0
            ranks = ";".join(read.classification[1]) if read.classification else ""
            outf_complete.write(
                f"{lineage}\t{read.n_reads}\t{read_frac}\t"
                f"{read.n_nucleotides}\t{ranks}\n"
            )


def make_bin_table(settings: RatSettings, rat_result: RatResult) -> None:
    total_reads = rat_result.mapping.sum_of_reads
    divisor = get_proper_fraction(rat_result.bin_data)

    with settings.files.bin_abundance.open("w", encoding="utf-8") as op:
        op.write(f"## command: {' '.join(sys.argv)}\n")
        op.write(
            "# bin\tnumber of reads\tfraction of reads\t"
            "bin length\tcorrected fraction\n"
        )

        for bin_id, read in rat_result.bin_data.items():
            read_fraction = read.n_reads / total_reads if total_reads else 0.0
            corrected_frac = (
                read.n_reads / read.n_nucleotides / divisor
                if read.n_nucleotides and divisor else 0.0
            )
            op.write(
                f"{bin_id}\t{read.n_reads}\t{read_fraction}\t"
                f"{read.n_nucleotides}\t{corrected_frac}\n"
            )


def write_read_table(settings: RatSettings, rat_result: RatResult) -> None:
    with settings.files.read_classification.open("w", encoding="utf-8") as outf:
        outf.write(f"## command: {' '.join(sys.argv)}\n")
        outf.write(
            "# read\tclassification\tbin classification\t"
            "contig classification\tdirect_mapping_contig\t"
            "direct_mapping_read\n"
            #'[aligner] unclassified contig\t[aligner] unclassified read\n'
        )

        for name, read in rat_result.mapping.read_dict.items():
            for c in range(len(read["contig"])):
                b_taxon = read["taxon_bin"][c]
                c_taxon = read["taxon_contig"][c]
                c_dm_taxon = (
                    "" if b_taxon or c_taxon
                    else read["taxon_contig_dm"][c]
                )
                if b_taxon is not None or c_taxon is not None or c_dm_taxon is not None:
                    assigned = "taxid assigned"
                else:
                    assigned = "no taxid assigned"
                outf.write(
                    f"{name}\t{assigned}\t{b_taxon}\t"
                    f"{c_taxon}\t{c_dm_taxon}\t\n"
                )

            for r_dm_taxon in read["taxon_read_dm"]:
                if r_dm_taxon is not None:
                    assigned = "taxid assigned"
                else:
                    assigned = "no taxid assigned"
                outf.write(f"{name}\t{assigned}\t\t\t\t{r_dm_taxon}\n")


def write_reports(settings: RatSettings, report: Report, rat_result: RatResult) -> RatResult:
    make_tax_table(settings, rat_result)
    report.progress_next()

    if settings.files.bin_abundance is not None:
        make_bin_table(settings, rat_result)
    report.progress_next()

    if settings.files.read_classification is not None:
        write_read_table(settings, rat_result)
    report.progress_next()
    return rat_result

    return
if __name__ == '__main__':
    sys.exit('Run \'CAT_pack reads\' to run Read Annotation Tool (RAT).')
