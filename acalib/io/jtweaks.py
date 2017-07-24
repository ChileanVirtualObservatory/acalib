from astropy.units.format.latex import Latex
import numpy as np
import astropy.visualization
from IPython import get_ipython


def jovial_array_styler(arr):
    lrep = getattr(arr, "_repr_latex_", None)
    if callable(lrep):
        return arr._repr_latex_()[1:-1]
    return repr(arr)


def jovial_array_makeup(a):
    lrep = getattr(a, "_repr_latex_", None)
    if callable(lrep):
        return a._repr_latex_()
    latex_value = np.array2string(
        a.view(),
        formatter=formatter,
        style=(Latex.format_exponential_notation
               if a.dtype.kind == 'f' else repr),
        max_line_width=np.inf,
        separator=',~')
    latex_value = latex_value.replace('...', r'\dots')
    return '${0}$'.format(latex_value)


def jovial_tuple_makeup(t):
    latex_value = ""
    for elem in t:
        if latex_value != "":
            latex_value += ",~"
        if isinstance(elem, tuple):
            latex_value += jovial_tuple_makeup(elem)[1:-1]
        elif isinstance(elem, np.ndarray):
            latex_value += jovial_array_makeup(elem)[1:-1]
        else:
            latex_value += Latex.format_exponential_notation(elem)
    return '$({0})$'.format(latex_value)

try:
    fmt = get_ipython().display_formatter.formatters['text/latex']
    astropy.visualization.quantity_support()
    formatter = {'float_kind': Latex.format_exponential_notation, 'numpystr': jovial_array_styler}
    np.set_printoptions(threshold=astropy.units.quantity.conf.latex_array_threshold)
    fmt.for_type(np.array, jovial_array_makeup)
    fmt.for_type(np.ndarray, jovial_array_makeup)
    #fmt.for_type(tuple, jovial_tuple_makeup)
except (NameError, AttributeError) as e :
    astropy.log.info("No IPython detected. Using standalone Python mode.")
    pass
