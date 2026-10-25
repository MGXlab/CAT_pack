""""configuration module of CAT pack

This module consist of three files `options.py`, `settings.py` and `validation.py`

Options.py
This python script holds dataclasses which can be filled with user input via the cli. Default values that the user can edit are also settable here.

Settings.py
This python script holds dataclasses which have been filled with validated input via validation.py. Default values that shouldn't be changed by the user via
the CLI can be declared here. As well as commands to run certain external tools. 

Validation.py
This script validates user input. Checks if the user needs to change values, files exist, files are empty or correct non corrupt input is supplies. 
and gives back guided errors with custom erros from `utils/errors.py' if input isn't as expected. This script depends on the python script `checks.py` to perfrom checks
"""
