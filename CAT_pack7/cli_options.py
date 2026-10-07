from functools import wraps
from inspect import Parameter, Signature, signature
from pathlib import Path
from typing import Annotated

from typer import Option

from .config.options import DiamondOptions, ExecutionOptions, MMseqsOptions

DiamondPathOption = Annotated[Path | None, Option(
    "--path-to-diamond", rich_help_panel="DIAMOND",
    help="Path to DIAMOND. Supply if it is not on PATH.",
)]


def execution_options(
    threads: Annotated[int, Option("--threads", "-n", min=1)] = ExecutionOptions.threads,
    verbose: Annotated[bool, Option("--verbose", help="Show aligner stdout.")] = ExecutionOptions.verbose,
    quiet: Annotated[bool, Option("--quiet", help="Turns off logging in the terminal.")] = ExecutionOptions.quiet,
    debug: Annotated[bool, Option("--debug", help="Show unexpected-error tracebacks.")] = ExecutionOptions.debug,
) -> ExecutionOptions:
    return ExecutionOptions(threads=threads, verbose=verbose, quiet=quiet, debug=debug)


def diamond_options(
    diamond_mode: Annotated[str, Option(
        "--diamond-mode", rich_help_panel="DIAMOND",
        help="default, faster, fast, mid-sensitive, sensitive, "
             "more-sensitive, very-sensitive, ultra-sensitive",
    )] = DiamondOptions.mode,
    block_size: Annotated[float, Option(
        "--block-size", rich_help_panel="DIAMOND",
        help="DIAMOND block-size. Lower uses less RAM/tmp.",
    )] = DiamondOptions.block_size,
    index_chunks: Annotated[int, Option(
        "--index-chunks", min=1, rich_help_panel="DIAMOND",
        help="Set to 4 on low-memory machines.",
    )] = DiamondOptions.index_chunks,
    no_self_hits: Annotated[bool, Option(
        "--no-self-hits", rich_help_panel="DIAMOND",
        help="Do not report identical self hits by DIAMOND.",
    )] = DiamondOptions.no_self_hits,
    path_to_diamond: DiamondPathOption = DiamondOptions.path_to_diamond,
) -> DiamondOptions:
    return DiamondOptions(
        mode=diamond_mode, block_size=block_size, index_chunks=index_chunks,
        no_self_hits=no_self_hits, path_to_diamond=path_to_diamond,
    )


def mmseqs_options(
    sensitivity: Annotated[float, Option(
        "--sensitivity", min=1.0, max=7.5, rich_help_panel="MMseqs2",
        help="MMseqs2 sensitivity (-s).",
    )] = MMseqsOptions.sensitivity,
    split_memory_limit: Annotated[str, Option(
        "--split-memory-limit", rich_help_panel="MMseqs2",
        help="MMseqs2 max memory per split, e.g. 10M, 1G. 0 uses all available memory.",
    )] = MMseqsOptions.split_memory_limit,
    path_to_mmseqs: Annotated[Path | None, Option(
        "--path-to-mmseqs", rich_help_panel="MMseqs2",
        help="Path to MMseqs2. Supply if it is not on PATH.",
    )] = MMseqsOptions.executable,
) -> MMseqsOptions:
    return MMseqsOptions(
        sensitivity=sensitivity, split_memory_limit=split_memory_limit,
        executable=path_to_mmseqs,
    )





def with_option_groups(**providers):
    """In-house dependency injection for CLI option groups.

    Use below @app.command()!!!
    @with_option_groups(execution=execution_options)

    You can provide multiple option groups
    """

    def decorate(command):
        command_signature = signature(command, eval_str=True)
        provider_signatures = {
            provider_name: signature(provider_func, eval_str=True)
            for provider_name, provider_func in providers.items()
        }

        command_parameters = [
            param for param_name, param in command_signature.parameters.items()
            if param_name not in providers
        ]
        provider_parameters = [
            param.replace(kind=Parameter.KEYWORD_ONLY)
            for provider_signature in provider_signatures.values()
            for param in provider_signature.parameters.values()
        ]
        merged_signature = command_signature.replace(parameters=[*command_parameters, *provider_parameters])

        @wraps(command)
        def wrapper(**kwargs):
            bound_arguments = merged_signature.bind(**kwargs)
            bound_arguments.apply_defaults()
            argument_values = dict(bound_arguments.arguments)

            injected_groups = {
                provider_name: providers[provider_name](**{
                    param_name: argument_values.pop(param_name)
                    for param_name in provider_signature.parameters
                })
                for provider_name, provider_signature in provider_signatures.items()
            }
            return command(**argument_values, **injected_groups)

        wrapper.__signature__ = merged_signature
        wrapper.__annotations__ = {
            param.name: param.annotation
            for param in merged_signature.parameters.values()
            if param.annotation is not Parameter.empty
        }
        if merged_signature.return_annotation is not Signature.empty:
            wrapper.__annotations__["return"] = merged_signature.return_annotation
        return wrapper

    return decorate
