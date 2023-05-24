#!/usr/bin/env python3

import argparse
import decimal
import multiprocessing
import os
import sys

import about
import check
import shared
import tax


def parse_arguments():
    parser = argparse.ArgumentParser(
            prog="CAT bins",
            description="Run Bin Annotation Tool (BAT) on a set of bins.",
            usage="CAT bins -b DIR -d DIR -t DIR [options] [-h / --help]",
            add_help=False)

    required = parser.add_argument_group("Required arguments")
    shared.add_argument(required, "bin_folder", True)
    shared.add_argument(required, "database_folder", True)
    shared.add_argument(required, "taxonomy_folder", True)

    optional = parser.add_argument_group("Optional arguments")
    shared.add_argument(optional, "bin_suffix", False, default=".fna")
    shared.add_argument(optional, "r", False, default=decimal.Decimal(5))
    shared.add_argument(optional, "f", False, default=decimal.Decimal(0.3))
    shared.add_argument(optional, "out_prefix", False, default="./out.BAT")
    shared.add_argument(optional, "proteins_fasta", False,
            help_=(
                "Path to concatenated predicted proteins fasta file "
                "generated during an earlier run of BAT on the same bins. If "
                "supplied, BAT will skip the protein prediction step."))
    shared.add_argument(optional, "alignment_file", False,
            help_=(
                "Path to alignment table generated during an earlier run of "
                "BAT on the same bins. If supplied, BAT will skip the "
                "alignment step and directly classify the bins. A "
                "concatenated predicted proteins fasta file should also be "
                "supplied with argument [-p / --proteins]."))
    shared.add_argument(optional, "path_to_prodigal", False,
            default="prodigal")
    shared.add_argument(optional, "path_to_diamond", False, default="diamond")
    shared.add_argument(optional, "no_stars", False)
    shared.add_argument(optional, "force", False)
    shared.add_argument(optional, "quiet", False)
    shared.add_argument(optional, "verbose", False)
    shared.add_argument(optional, "no_log", False)
    shared.add_argument(optional, "help", False)
    shared.add_argument(optional, "IkwId", False)

    specific = parser.add_argument_group("DIAMOND specific optional arguments")
    shared.add_all_diamond_arguments(specific)

    (args, extra_args) = parser.parse_known_args()

    extra_args = [arg for (i, arg) in enumerate(extra_args) if
                  (i, arg) != (0, "bins")]
    if len(extra_args) > 0:
        sys.exit("error: too much arguments supplied:\n{0}".format(
            "\n".join(extra_args)))
        
    # Check experimental features.
    if not args.IkwId:
        if args.top < 15:
            sys.exit("error: [--top] can only be set lower than 15 with the "
                    "[--I_know_what_Im_doing] flag. See README.md as to why "
                    "this is the case.")
            
        if args.r > 15 and args.alignment_file:
            sys.exit("error: [-r / --range] can only be set higher than 15 in "
                    "combination with [-a / --diamond_alignment] with the "
                    "[--I_know_what_Im_doing] flag. See README.md as to why "
                    "this is the case.")
            
    # Add extra arguments.
    shared.expand_arguments(args)
            
    return args


def import_bins(bin_folder, bin_suffix, log_file, quiet):
    message = "Importing bins from {0}.".format(bin_folder)
    shared.give_user_feedback(message, log_file, quiet)

    bin2contigs = {}
    contig_names = set()

    for file_ in os.listdir(bin_folder):
        if file_.startswith("."):
            # Skip hidden files.
            continue

        if not file_.endswith(bin_suffix):
            continue
        
        if ".concatenated." in file_:
            # Skip concatenated contig fasta and predicted protein fasta files
            # from earlier runs.
            continue
        
        # Keep the suffix in the bin name.
        bin_ = file_

        bin2contigs[bin_] = []

        with open("{0}{1}".format(bin_folder, file_), "r") as f1:
            for line in f1:
                if line.startswith(">"):
                    contig = line.split()[0].rstrip().lstrip(">")

                    # Add bin name in front of the contig name.
                    new_contig_name = "{0}_{1}".format(bin_, contig)
                    
                    if new_contig_name in contig_names:
                        message = (
                                "BAT has encountered {0} twice in bin {1}. "
                                "Each fasta header should be unique in each "
                                "bin.".format(contig, bin_))
                        shared.give_user_feedback(
                                message, log_file, quiet, error=True)

                        sys.exit(1)

                    contig_names.add(new_contig_name)

                    bin2contigs[bin_].append(new_contig_name)
                    
    if len(bin2contigs) == 1:
        message = "1 bin found!"
    else:
        message = "{0:,d} bins found!".format(len(bin2contigs))
    shared.give_user_feedback(message, log_file, quiet)

    return (bin2contigs, contig_names)


def make_concatenated_fasta(
        concatenated_fasta, bin2contigs, bin_folder, log_file, quiet):
    message = "Writing {0}.".format(concatenated_fasta)
    shared.give_user_feedback(message, log_file, quiet)

    with open(concatenated_fasta, "w") as outf1:
        for bin_ in sorted(bin2contigs):
            with open("{0}{1}".format(bin_folder, bin_), "r") as f1:
                for line in f1:
                    if line.startswith(">"):
                        contig = line.split()[0].rstrip().lstrip(">")
                        
                        # add bin name in front of the contig name.
                        outf1.write(">{0}_{1}\n".format(bin_, contig))
                    else:
                        outf1.write(line)

    return
                        
                        
def run():
    args = parse_arguments()

    message = "# CAT v{0}.".format(about.__version__)
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)
    
    # Check at which state to start.
    step_list = []
    if not args.proteins_fasta and not args.alignment_file:
        message = (
                "\n"
                "BAT is running. Protein prediction, alignment, and bin "
                "classification are carried out.")
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=False)

        step_list.append("predict_proteins")
        step_list.append("align")
    elif args.proteins_fasta and not args.alignment_file:
        message = (
                "\n"
                "BAT is running. Since a predicted protein fasta is supplied, "
                "only alignment and bin classification are carried out.")
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=False)

        step_list.append("align")
    elif args.proteins_fasta and args.alignment_file:
        message = (
                "\n"
                "BAT is running. Since a predicted protein fasta and "
                "alignment file are supplied, only bin classification is "
                "carried out.")
        shared.give_user_feedback(message, args.log_file, args.quiet,
                show_time=False)
    elif not args.proteins_fasta and args.alignment_file:
        message = (
                "if you want BAT to directly classify a set of bins, you "
                "should not only supply a DIAMOND alignment table but also a "
                "concatenated predicted protein fasta file with argument "
                "[-p / --proteins].")
        shared.give_user_feedback(message, args.log_file, args.quiet,
                error=True)

        sys.exit(1)

    step_list.append("classify")

    # Print variables.
    message = (
            "Rarw!\n\n"
            "Supplied command: {0}\n\n"
            "Bin folder: {1}\n"
            "Taxonomy folder: {2}\n"
            "Database folder: {3}\n"
            "Parameter r: {4}\n"
            "Parameter f: {5}\n"
            "Log file: {6}\n\n"
            "-----------------\n".format(
                " ".join(sys.argv),
                args.bin_folder,
                args.taxonomy_folder,
                args.database_folder,
                int(args.r),
                float(args.f),
                args.log_file))
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)

    # Check binaries, output files, taxonomy folder and database folder, and
    # set variables.
    message = "Doing some pre-flight checks first."
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)

    errors = []

    errors.append(
            check.check_bin_folder(
                args.bin_folder, args.bin_suffix, args.log_file, args.quiet))
    
    errors.append(
            check.check_out_prefix(args.out_prefix, args.log_file, args.quiet))
    
    if "predict_proteins" in step_list:
        errors.append(
                check.check_prodigal_binaries(
                    args.path_to_prodigal, args.log_file, args.quiet))

        setattr(args,
                "concatenated_fasta",
                "{0}.concatenated.fasta".format(args.out_prefix))
        setattr(args,
                "proteins_fasta",
                "{0}.concatenated.predicted_proteins.faa".format(
                    args.out_prefix))
        setattr(args,
                "proteins_gff",
                "{0}.concatenated.predicted_proteins.gff".format(
                    args.out_prefix))

        if not args.force:
            errors.append(
                    check.check_output_file(
                        args.concatenated_fasta, args.log_file, args.quiet))
            errors.append(
                    check.check_output_file(
                        args.proteins_fasta, args.log_file, args.quiet))
            errors.append(
                    check.check_output_file(
                        args.proteins_gff, args.log_file, args.quiet))
            
    if "align" in step_list:
        errors.append(
                check.check_diamond_binaries(
                    args.path_to_diamond, args.log_file, args.quiet))

        setattr(args,
                "alignment_file",
                "{0}.concatenated.alignment.diamond".format(args.out_prefix))

        if not args.force:
            errors.append(
                    check.check_output_file(
                        args.alignment_file, args.log_file, args.quiet))

    errors.append(
            check.check_folders_for_run(
                args.taxonomy_folder,
                args.nodes_dmp,
                args.names_dmp,
                args.database_folder,
                args.diamond_database,
                args.fastaid2LCAtaxid_file,
                args.taxids_with_multiple_offspring_file,
                step_list,
                args.log_file,
                args.quiet))

    setattr(args,
            "bin2classification_output_file",
            "{0}.bin2classification.txt".format(args.out_prefix))
    setattr(args,
            "ORF2LCA_output_file",
            "{0}.ORF2LCA.txt".format(args.out_prefix))

    if not args.force:
        errors.append(
                check.check_output_file(
                    args.bin2classification_output_file,
                    args.log_file,
                    args.quiet))
        errors.append(
                check.check_output_file(
                    args.ORF2LCA_output_file, args.log_file, args.quiet))
        
    if "predict_proteins" not in step_list:
        errors.append(
                check.check_fasta(
                    args.proteins_fasta, args.log_file, args.quiet))

    if "align" in step_list:
        errors.append(
                check.check_top(args.top, args.r, args.log_file, args.quiet))

    # Print all variables.
    shared.print_variables(args, step_list)

    if True in errors:
        sys.exit(1)

    message = "Ready to fly!\n\n-----------------\n"
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)
    
    # Start BAT.
    (bin2contigs, contig_names) = import_bins(
            args.bin_folder, args.bin_suffix, args.log_file, args.quiet)

    if "predict_proteins" in step_list:
        make_concatenated_fasta(
                args.concatenated_fasta,
                bin2contigs,
                args.bin_folder,
                args.log_file,
                args.quiet)

        shared.run_prodigal(
                args.path_to_prodigal,
                args.concatenated_fasta,
                args.proteins_fasta,
                args.proteins_gff,
                args.log_file,
                args.quiet)
        
    contig2ORFs = shared.import_ORFs(
            args.proteins_fasta, args.log_file, args.quiet)
    
    check.check_whether_ORFs_are_based_on_contigs(
            contig_names, contig2ORFs, args.log_file, args.quiet)
    
    if "align" in step_list:
        shared.run_diamond(args)

    (ORF2hits,
            all_hits) = shared.parse_tabular_alignment(
                    args.alignment_file,
                    args.one_minus_r,
                    args.log_file,
                    args.quiet)

    (taxid2parent,
            taxid2rank) = tax.import_nodes(
            args.nodes_dmp, args.log_file, args.quiet)
    fastaid2LCAtaxid = tax.import_fastaid2LCAtaxid(
            args.fastaid2LCAtaxid_file, all_hits, args.log_file, args.quiet)
    taxids_with_multiple_offspring = tax.import_taxids_with_multiple_offspring(
            args.taxids_with_multiple_offspring_file,
            args.log_file,
            args.quiet)
    
    message = "BAT is flying! Files {0} and {1} are created.".format(
        args.bin2classification_output_file, args.ORF2LCA_output_file)
    shared.give_user_feedback(message, args.log_file, args.quiet)

    n_classified_bins = 0

    with open(args.bin2classification_output_file, "w") as outf1, open(args.ORF2LCA_output_file, "w") as outf2:
        outf1.write("# bin\tclassification\treason\tlineage\tlineage scores\n")

        outf2.write("# ORF\tbin\tnumber of hits\tlineage\ttop bit-score\n")
        
        for bin_ in sorted(bin2contigs):
            LCAs_ORFs = []

            for contig in sorted(bin2contigs[bin_]):
                if contig not in contig2ORFs:
                    continue

                for ORF in contig2ORFs[contig]:
                    if ORF not in ORF2hits:
                        outf2.write("{0}\t{1}\tORF has no hit to database\n"
                                "".format(ORF, bin_))

                        continue

                    n_hits = len(ORF2hits[ORF])

                    (taxid,
                            top_bitscore) = tax.find_LCA_for_ORF(
                                    ORF2hits[ORF],
                                    fastaid2LCAtaxid,
                                    taxid2parent)
                     
                    if taxid.startswith("no taxid found"):
                        outf2.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                            ORF, bin_, n_hits, taxid, top_bitscore))
                    else:
                        lineage = tax.find_lineage(taxid, taxid2parent)

                        if not args.no_stars:
                            lineage = tax.star_lineage(
                                    lineage, taxids_with_multiple_offspring)

                        outf2.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                            ORF,
                            bin_,
                            n_hits,
                            ";".join(lineage[::-1]),
                            top_bitscore))
                                       
                    LCAs_ORFs.append((taxid, top_bitscore),)
                    
            if len(LCAs_ORFs) == 0:
                outf1.write("{0}\tno taxid assigned\tno hits to database\n"
                        "".format(bin_))

                continue

            (lineages,
                    lineages_scores,
                    based_on_n_ORFs) = tax.find_weighted_LCA(
                            LCAs_ORFs, taxid2parent, args.f)

            if lineages == "no ORFs with taxids found.":
                outf1.write("{0}\tno taxid assigned\t"
                        "hits not found in taxonomy files\n".format(bin_))

                continue

            if lineages == "no lineage whitelisted.":
                outf1.write(
                        "{0}\tno taxid assigned\t"
                        "no lineage reached minimum bit-score support\n"
                        "".format(bin_))

                continue
            
            # The bin has a valid classification.
            n_classified_bins += 1

            total_n_ORFs = sum([len(contig2ORFs[contig]) for
                contig in bin2contigs[bin_] if contig in contig2ORFs])
            
            for (i, lineage) in enumerate(lineages):
                if not args.no_stars:
                    lineage = tax.star_lineage(
                            lineage, taxids_with_multiple_offspring)
                
                scores = ["{0:.2f}".format(score) for
                        score in lineages_scores[i]]
                
                if len(lineages) == 1:
                    # There is only one classification.
                    outf1.write(
                            "{0}\t"
                            "taxid assigned\t"
                            "based on {1}/{2} ORFs\t"
                            "{3}\t"
                            "{4}\n".format(
                                bin_,
                                based_on_n_ORFs,
                                total_n_ORFs,
                                ";".join(lineage[::-1]),
                                ";".join(scores[::-1])))
                else:
                    # There are multiple classifications.
                    outf1.write(
                            "{0}\t"
                            "taxid assigned ({1}/{2})\t"
                            "based on {3}/{4} ORFs\t"
                            "{5}\t"
                            "{6}\n".format(
                                bin_,
                                i + 1,
                                len(lineages),
                                based_on_n_ORFs,
                                total_n_ORFs,
                                ";".join(lineage[::-1]),
                                ";".join(scores[::-1])))
                                   
    message = ("\n-----------------\n\n"
            "{0} BAT is done! {1:,d}/{2:,d} bins have taxonomy assigned."
            "".format(shared.timestamp(), n_classified_bins, len(bin2contigs)))
    shared.give_user_feedback(message, args.log_file, args.quiet,
            show_time=False)
  
    if args.f < 0.5:
        message = ("since f is set to smaller than 0.5, one bin "
                "may have multiple classifications.")
        shared.give_user_feedback(
                message,
                args.log_file,
                args.quiet,
                show_time=False,
                warning=True
                )

    return


if __name__ == "__main__":
    sys.exit("Run \'CAT bins\' to run Bin Annotation Tool (BAT) on a "
             "set of bins.")
