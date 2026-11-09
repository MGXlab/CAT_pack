"""External tool management module of CAT pack

This module consist of four files `aligner.py`. `bwa_mem.py`, `pyrodigal.py`, `tool_runner.py`

aligner.py
This python script runs the protein aligner DIAMOND, or MMseqs2.

gpu.py
This python script checks for CUDA and MMseqs2 compatibility.

mapper.py
This python script runs creates a sorted bam file. It uses bwa mem or bwa-mem2
and pipes the output into samtools directly

pyrodigal.py
runs protein prediction (with now only pyrodigal being supported).

tool_runner.py
The tool runner, this handles the drainage of log output for external tools.
In terminal (visible with --debug) and for the .log file

"""
