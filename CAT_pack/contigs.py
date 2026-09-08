#!/usr/bin/env python3

import argparse
import decimal
import multiprocessing
import sys

import about
import check
import shared
import tax
import classification


def parse_arguments():
    parser = argparse.ArgumentParser(
            prog="CAT_pack contigs",
            description="Run Contig Annotation Tool (CAT).",
            usage=("CAT_pack contigs -c <FILE> -d <DIR> -t <DIR> [options] "
                "[-h / --help]"),
            add_help=False
            )
    
    required = parser.add_argument_group("Required arguments")
    shared.add_argument(required, "contigs_fasta", True)
    shared.add_argument(required, "database_folder", True)
    shared.add_argument(required, "taxonomy_folder", True)

    optional = parser.add_argument_group("Optional arguments")
    shared.add_argument(optional, "r", False, default=decimal.Decimal(10))
    shared.add_argument(optional, "f", False, default=decimal.Decimal(0.5))
    shared.add_argument(optional, "out_prefix", False, default="./out.CAT")
    shared.add_argument(optional, "proteins_fasta", False)
    shared.add_argument(optional, "alignment_file", False)
    shared.add_argument(optional, "aligner", False, default="DIAMOND")
    shared.add_argument(optional, "no_stars", False)
    shared.add_argument(optional, "top", False, default=11)
    shared.add_argument(
            optional, "nproc", False, default=multiprocessing.cpu_count())
    shared.add_argument(optional, "compress", False)
    shared.add_argument(optional, "tmpdir", False)
    shared.add_argument(optional, "force", False)
    shared.add_argument(optional, "quiet", False)
    shared.add_argument(optional, "verbose", False)
    shared.add_argument(optional, "no_log", False)
    shared.add_argument(optional, "help", False)
    shared.add_argument(optional, "IkwId", False)

    specific = parser.add_argument_group("DIAMOND specific optional arguments")
    shared.add_all_diamond_arguments(specific)

    specific = parser.add_argument_group("MMseqs2 specific optional arguments")
    shared.add_all_mmseqs2_arguments(specific)

    args, extra_args = parser.parse_known_args()
    
    extra_args = [arg for (i, arg) in enumerate(extra_args) if
            (i, arg) != (0, "contigs")]
    if len(extra_args) > 0:
        sys.exit("error: too many arguments supplied:\n{0}".format(
            "\n".join(extra_args)))
        
    # Check experimental features.
    if not args.IkwId:
        if args.top < 11:
            sys.exit(
                    "error: --top can only be set lower than 11 with the "
                    "--I_know_what_Im_doing flag. See README.md as to why "
                    "this is the case."
                    )
            
        if args.r > 11 and args.alignment_file:
            sys.exit(
                    "error: --range can only be set higher than 11 in "
                    "combination with --alignment_table with the "
                    "--I_know_what_Im_doing flag. See README.md as to why "
                    "this is the case."
                    )

    # Add extra arguments.
    shared.expand_arguments(args)

    return args


def run():
    args = parse_arguments()

    message = "# CAT_pack v{0}.".format(about.__version__)
    shared.give_user_feedback(
            message, args.log_file, args.quiet, show_time=False)

    # Check at which state to start.
    step_list = []
    if not args.proteins_fasta and not args.alignment_file:
        message = (
                "\n"
                "CAT is running. Protein prediction, alignment, and contig "
                "classification are carried out."
                )
        shared.give_user_feedback(
                message, args.log_file, args.quiet, show_time=False)

        step_list.append("predict_proteins")
        step_list.append("align")
    elif args.proteins_fasta and not args.alignment_file:
        message = (
                "\n"
                "CAT is running. Since a predicted protein fasta is supplied, "
                "only alignment and contig classification are carried out."
                )
        shared.give_user_feedback(
                message, args.log_file, args.quiet, show_time=False)

        step_list.append("align")
    elif args.proteins_fasta and args.alignment_file:
        message = (
                "\n"
                "CAT is running. Since a predicted protein fasta and "
                "alignment file are supplied, only contig classification is "
                "carried out."
                )
        shared.give_user_feedback(
                message, args.log_file, args.quiet, show_time=False)
    elif not args.proteins_fasta and args.alignment_file:
        message = (
                "if you want CAT to directly do the classification, you "
                "should not only supply an alignment table but also a "
                "predicted protein fasta file with argument --proteins_fasta."
                )
        shared.give_user_feedback(
                message, args.log_file, args.quiet, error=True)

        sys.exit(1)

    step_list.append("classify")

    # Print variables.
    message = (
            "Rarw!\n\n"
            "Supplied command: {0}\n\n"
            "Contigs fasta: {1}\n"
            "Taxonomy folder: {2}\n"
            "Database folder: {3}\n"
            "Parameter r: {4}\n"
            "Parameter f: {5}\n"
            "Log file: {6}\n\n"
            "-----------------\n".format(
                " ".join(sys.argv),
                args.contigs_fasta,
                args.taxonomy_folder,
                args.database_folder,
                int(args.r),
                float(args.f),
                args.log_file)
            )
    shared.give_user_feedback(
            message, args.log_file, args.quiet, show_time=False)

    # Check binaries, output files, taxonomy folder and database folder, and
    # set variables.
    message = "Doing some pre-flight checks first."
    shared.give_user_feedback(
            message, args.log_file, args.quiet, show_time=False)

    errors = []

    errors.append(
            check.check_out_prefix(args.out_prefix, args.log_file, args.quiet))
    
    if "predict_proteins" in step_list:
        errors.append(
                check.check_pyrodigal_install(
                    args.log_file, args.quiet)
                )

        setattr(
                args,
                "proteins_fasta",
                "{0}.predicted_proteins.faa".format(args.out_prefix)
                )
        setattr(
                args,
                "proteins_gff",
                "{0}.predicted_proteins.gff".format(args.out_prefix)
                )

        if not args.force:
            errors.append(
                    check.check_output_file(
                        args.proteins_fasta, args.log_file, args.quiet)
                    )
            errors.append(
                    check.check_output_file(
                        args.proteins_gff, args.log_file, args.quiet)
                    )
            
    if "align" in step_list:
        if args.aligner.lower() == "diamond":
            errors.append(
                    check.check_diamond_binaries(
                        args.path_to_diamond, args.log_file, args.quiet)
                    )

            setattr(
                    args,
                    "alignment_file",
                    "{0}.alignment.diamond".format(args.out_prefix)
                    )
        elif args.aligner.lower() == "mmseqs2":
            errors.append(
                    check.check_mmseqs2_binaries(
                        args.path_to_mmseqs2, args.log_file, args.quiet)
                    )
            setattr(
                    args,
                    "alignment_file",
                    "{0}.alignment.mmseqs2".format(args.out_prefix)
                    )
        else:
            # For debugging...
            sys.exit("Something wrong!")

        if not args.force:
            errors.append(
                    check.check_output_file(
                        args.alignment_file, args.log_file, args.quiet)
                    )

    errors.append(
            check.check_folders_for_run(
                args.taxonomy_folder,
                args.nodes_dmp,
                args.names_dmp,
                args.database_folder,
                args.aligner,
                args.diamond_database,
                args.mmseqs2_database,
                args.fastaid2LCAtaxid_file,
                args.taxids_with_multiple_offspring_file,
                step_list,
                args.log_file,
                args.quiet
                )
            )

    setattr(
            args,
            "contig2classification_output_file",
            "{0}.contig2classification.txt".format(args.out_prefix)
            )
    setattr(
            args,
            "ORF2LCA_output_file",
            "{0}.ORF2LCA.txt".format(args.out_prefix)
            )

    if not args.force:
        errors.append(
                check.check_output_file(
                    args.contig2classification_output_file,
                    args.log_file,
                    args.quiet
                    )
                )
        errors.append(
                check.check_output_file(
                    args.ORF2LCA_output_file, args.log_file, args.quiet)
                )

    if "predict_proteins" not in step_list:
        errors.append(
                check.check_fasta(
                    args.proteins_fasta, args.log_file, args.quiet)
                )

    if "align" in step_list:
        if not args.force:
            errors.append(
                    check.check_output_file(
                        args.alignment_file,
                        args.log_file,
                        args.quiet
                        )
                    )

        errors.append(
                check.check_top(args.top, args.r, args.log_file, args.quiet))

    # Print all variables.
    shared.print_variables(args, step_list)

    if True in errors:
        sys.exit(1)

    message = "Ready to fly!\n\n-----------------\n"
    shared.give_user_feedback(
            message, args.log_file, args.quiet, show_time=False)
    
    # Start CAT.
    contig_names = shared.import_contig_names(
            args.contigs_fasta, args.log_file, args.quiet)
    
    if "predict_proteins" in step_list:
        shared.run_pyrodigal(
                args.contigs_fasta,
                args.proteins_fasta,
                args.proteins_gff,
                args.nproc,
                args.log_file,
                args.quiet
                )
        
    contig2ORFs = shared.import_ORFs(
            args.proteins_fasta, args.log_file, args.quiet)
    
    check.check_whether_ORFs_are_based_on_contigs(
            contig_names, contig2ORFs, args.log_file, args.quiet)
    
    if "align" in step_list:
        shared.run_aligner(args)

    ORF2hits, all_hits = shared.parse_tabular_alignment(
            args.alignment_file, args.one_minus_r, args.log_file, args.quiet)

    taxid2parent, taxid2rank = tax.import_nodes(
            args.nodes_dmp, args.log_file, args.quiet)
    fastaid2LCAtaxid = tax.import_fastaid2LCAtaxid(
            args.fastaid2LCAtaxid_file, all_hits, args.log_file, args.quiet)
    taxids_with_multiple_offspring = tax.import_taxids_with_multiple_offspring(
            args.taxids_with_multiple_offspring_file,
            args.log_file,
            args.quiet
            )



    message = "CAT is spinning! Files {0} and {1} are created.".format(
            args.contig2classification_output_file, args.ORF2LCA_output_file)
    shared.give_user_feedback(message, args.log_file, args.quiet)

    cat_engine = classification.ClassificationEngine(
        taxid2parent=taxid2parent,
        fastaid2taxid=fastaid2LCAtaxid,
        fraction=args.f,
    )

    n_classified_contigs = 0
    
    with (
            open(args.contig2classification_output_file, "w") as outf1,
            open(args.ORF2LCA_output_file, "w") as outf2
            ):
        outf1.write(f"# contig\tclassification\treason\tlineage\t"
                f"lineage scores (f: {float(args.f)})\n")

        outf2.write(f"# ORF\tnumber of hits (r: {args.r})\tlineage\ttop bit-score\n")
        
        for contig in sorted(contig_names):

            result = cat_engine.classify_group(entity_id=contig,orf_ids=contig2ORFs.get(contig, ()), orf2hits=ORF2hits)

            for orf_result in result.orf_results:
                if orf_result.status == classification.ORFStatus.NO_HIT:
                    outf2.write(f"{orf_result.orf_id}\tORF has no hit to database\n")

                    continue

                if orf_result.status == classification.ORFStatus.NO_TAXID:
                    outf2.write(
                        f"{orf_result.orf_id}\t{orf_result.n_hits}\t"
                        f"{orf_result.message}\t{orf_result.top_bitscore}\n"
                    )
                    continue

                lineage = list(orf_result.lineage)

                if not args.no_stars:
                    lineage = tax.star_lineage(
                        lineage, taxids_with_multiple_offspring)
                    
                outf2.write(f"{orf_result.orf_id}\t{orf_result.n_hits}\t"
                            f"{';'.join(lineage[::-1])}\t{orf_result.top_bitscore}\n")

            if result.status == classification.ClassificationStatus.NO_ORFS:
                outf1.write(f"{contig}\tno taxid assigned\tno ORFs found\n")
                continue

            if result.status == classification.ClassificationStatus.NO_HITS:
                outf1.write(f"{contig}\tno taxid assigned\tno hits to database\n")
                continue

            if result.status == classification.ClassificationStatus.NO_TAXIDS:
                outf1.write(f"{contig}\tno taxid assigned\thits not found in taxonomy files\n")
                continue

            if result.status == classification.ClassificationStatus.NO_LINEAGE_SUPPORT:
                outf1.write(f"{contig}\tno taxid assigned\t"
                            f"no lineage reached minimum bit-score support\n")
                continue

            # The contig has a valid classification.
            n_classified_contigs += 1

            for i, assignment in enumerate(result.assignments):
                lineage = list(assignment.lineage)
                if not args.no_stars:
                    lineage = tax.star_lineage(lineage,taxids_with_multiple_offspring)

                scores = [
                    f"{score:.2f}"
                    for score in assignment.lineage_scores
                ]

                if len(result.assignments) == 1:
                    outf1.write(
                        f"{contig}\t"
                        f"taxid assigned\t"
                        f"based on {result.based_on_n_ORFs}/{result.total_n_ORFs} ORFs\t"
                        f"{';'.join(lineage[::-1])}\t"
                        f"{';'.join(scores[::-1])}\n"
                    )

                else:
                    outf1.write(
                        f"{contig}\t"
                        f"taxid assigned ({i + 1}/{len(result.assignments)})\t"
                        f"based on {result.based_on_n_ORFs}/{result.total_n_ORFs} ORFs\t"
                        f"{';'.join(lineage[::-1])}\t"
                        f"{';'.join(scores[::-1])}\n"
                    )

    message = (
            "\n-----------------\n\n"
            "{0} CAT is done! {1:,d}/{2:,d} contigs ({3:.2f}%) have "
            "taxonomy assigned.".format(
                shared.timestamp(),
                n_classified_contigs,
                len(contig_names),
                n_classified_contigs / len(contig_names) * 100
                )
            )
    shared.give_user_feedback(
            message, args.log_file, args.quiet, show_time=False)

    if args.f < 0.5:
        message = ("since f is set to smaller than 0.5, one contig may have "
                "multiple classifications.")
        shared.give_user_feedback(
                message,
                args.log_file,
                args.quiet,
                show_time=False,
                warning=True
                )

    return


if __name__ == "__main__":
    sys.exit("Run \'CAT_pack contigs\' to run Contig Annotation Tool (CAT).")
