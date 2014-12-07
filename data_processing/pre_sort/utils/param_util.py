"""
Handle user-specified and default parameters from KK .prm file.
Copied from KlustaViewa 0.3.0.b1
"""
# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------
import os
import pprint

from six import string_types, iteritems




# -----------------------------------------------------------------------------
# Python script <==> dictionaries conversion
# -----------------------------------------------------------------------------
def get_params(filename=None, **kwargs):
    """Return all the parameters, retrieved following this order of priority:
    
    * parameters specified as keyword arguments in this function,
    * parameters specified in the .PRM file given in `filename`,
    * default parameters.
    
    """
    # Extract sample_rate before loading the default parameters.
    # This is because some default parameters are expressed as a function
    # of the sample rate.
    sample_rate = get_pydict(filename).get('sample_rate', None) or kwargs['sample_rate']
    if 'sample_rate' not in kwargs:
        kwargs['sample_rate'] = sample_rate
    default = load_default_params(kwargs)
    return get_pydict(filename=filename,
                      pydict_default=default,
                      **kwargs)


# -----------------------------------------------------------------------------
# Default parameters
# -----------------------------------------------------------------------------

def load_default_params(namespace=None):
    """Load default parameters, in a given namespace (empty by default)."""

    if namespace is None:
        namespace = {}
    # The default parameters are read in a namespace that must contain
    # sample_rate.
    assert namespace['sample_rate'] > 0
    folder = os.path.dirname(os.path.realpath(__file__))
    params_default_path = os.path.join(folder, 'params_default.py')
    with open(params_default_path, 'r') as f:
        params_default_python = f.read()
    params_default = python_to_pydict(params_default_python, namespace.copy())
    return to_lower(params_default)


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------
def display_params(prm):
    return pprint.pformat(prm)


def to_str(val):
    """Get a string representation of any Python variable."""
    if isinstance(val, string_types):
        return "'{0:s}'".format(val)
    else:
        return str(val)


def to_lower(d):
    return {key.lower(): val for key, val in iteritems(d)}


# -----------------------------------------------------------------------------
# Python script <==> dictionaries conversion
# -----------------------------------------------------------------------------
def python_to_pydict(script_contents, namespace=None):
    """Load a Python script with dictionaries into a dictionary."""
    if namespace is None:
        namespace = {}
    exec script_contents in {}, namespace
    return to_lower(namespace)


def pydict_to_python(pydict):
    """Convert a dictionaries dictionary into a Python script."""
    return "\n".join(["{0:s} = {1:s}".format(key, to_str(val))
                      for key, val in sorted(pydict.iteritems())])


def get_pydict(filename=None, pydict_default={}, **kwargs):
    """Load a Python script with key=value lines, convert it into a dictionary,
    and fill unset values to default values."""
    pydict_final = pydict_default.copy()
    if isinstance(filename, string_types):
        # Path to pydict file.
        with open(filename, 'r') as f:
            pydict_prm = python_to_pydict(f.read())
            pydict_final.update(pydict_prm)
    pydict_final.update(kwargs)
    return to_lower(pydict_final)
