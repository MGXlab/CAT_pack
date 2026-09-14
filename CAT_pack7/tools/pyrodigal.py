def run_protein_prediction(args, files, log, report, predictor):
    if predictor == "pyrodigal":
        return run_pyrodigal(args, files, log, report)
    return None


def run_pyrodigal(args, files, log, report):
    """placeholder for the pyrodigal"""
    pass