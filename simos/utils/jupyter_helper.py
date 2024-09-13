import numpy as np
import IPython.display
from ..trivial import flatten

def npdisp(a,atol=1e-10):
        """
        Generate a LaTeX representation a numpy array. Can be used for
        formatted output in a ipython notebook.
        Values smaller than atol will be treated as zero.
        
        """
        def _format_float(value,cutoff):
            #print("value", value, "cutoff", cutoff)
            if value == 0.0:
                return "0.0"
            elif abs(value) < cutoff:
                return "0.0"
            elif abs(value) > 1000.0 or abs(value) < 0.001:
                return ("%.3e" % value).replace("e", r"\times10^{") + "}"
            elif abs(value - int(value)) < 0.001:
                return "%.1f" % value
            else:
                return "%.3f" % value

        def _format_element(m, n, d):
            s = " & " if n > 0 else ""
            if isinstance(d, str):
                return s + d
            else:
                if abs(np.imag(d)) < atol:
                    return s + _format_float(np.real(d),atol)
                elif abs(np.real(d)) < atol:
                    return s + _format_float(np.imag(d),atol) + "i"
                else:
                    s_re = _format_float(np.real(d),atol)
                    s_im = _format_float(np.imag(d),atol)
                    if np.imag(d) > 0.0:
                        return (s + "(" + s_re + "+" + s_im + "j)")
                    else:
                        return (s + "(" + s_re + s_im + "j)")

        shape = np.shape(a)
        if len(shape) == 1:
            a = np.expand_dims(a,1)
            shape = np.shape(a)
        elif len(shape) > 2:
            raise ValueError('Only supporting 1D & 2D numpy arrays')
        s = r''

        s += r'\begin{equation*}\left(\begin{array}{*{11}c}'
        M, N = shape

        if M > 10 and N > 10:
            # truncated matrix output
            for m in range(5):
                for n in range(5):
                    s += _format_element(m, n, a[m, n])
                s += r' & \cdots'
                for n in range(N - 5, N):
                    s += _format_element(m, n, a[m, n])
                s += r'\\'

            for n in range(5):
                s += _format_element(m, n, r'\vdots')
            s += r' & \ddots'
            for n in range(N - 5, N):
                s += _format_element(m, n, r'\vdots')
            s += r'\\'

            for m in range(M - 5, M):
                for n in range(5):
                    s += _format_element(m, n, a[m, n])
                s += r' & \cdots'
                for n in range(N - 5, N):
                    s += _format_element(m, n, a[m, n])
                s += r'\\'

        elif M > 10 and N <= 10:
            # truncated vertically elongated matrix output
            for m in range(5):
                for n in range(N):
                    s += _format_element(m, n, a[m, n])
                s += r'\\'

            for n in range(N):
                s += _format_element(m, n, r'\vdots')
            s += r'\\'

            for m in range(M - 5, M):
                for n in range(N):
                    s += _format_element(m, n, a[m, n])
                s += r'\\'

        elif M <= 10 and N > 10:
            # truncated horizontally elongated matrix output
            for m in range(M):
                for n in range(5):
                    s += _format_element(m, n, a[m, n])
                s += r' & \cdots'
                for n in range(N - 5, N):
                    s += _format_element(m, n, a[m, n])
                s += r'\\'

        else:
            # full output
            for m in range(M):
                for n in range(N):
                    s += _format_element(m, n, a[m, n])
                s += r'\\'

        s += r'\end{array}\right)\end{equation*}'
        
        return IPython.display.Latex(s)

def bar(sequence, every=None, size=None, name='Items'):
    """Progress bar for jupyter notebooks
    
    Parameters:
    sequence : sequence to be stepped through

    Optional parameters:
    every    : update how often?
    name     : Naming in the display
    """

    from ipywidgets import IntProgress, HTML, VBox
    from IPython.display import display
    is_iterator = False
    if size is None:
        try:
            size = len(sequence)
        except TypeError:
            is_iterator = True
    if size is not None:
        if every is None:
            if size <= 200:
                every = 1
            else:
                every = int(size / 200)     # every 0.5%
    else:
        assert every is not None, 'sequence is iterator, set every'

    if is_iterator:
        progress = IntProgress(min=0, max=1, value=1)
        progress.bar_style = 'info'
    else:
        progress = IntProgress(min=0, max=size, value=0)
    label = HTML()
    box = VBox(children=[label, progress])
    display(box)

    index = 0
    try:
        for index, record in enumerate(sequence, 1):
            if index == 1 or index % every == 0:
                if is_iterator:
                    label.value = '{name}: {index} / ?'.format(
                        name=name,
                        index=index
                    )
                else:
                    progress.value = index
                    label.value = u'{name}: {index} / {size}'.format(
                        name=name,
                        index=index,
                        size=size
                    )
            yield record
    except:
        progress.bar_style = 'danger'
        raise
    else:
        progress.bar_style = 'success'
        progress.value = index
        label.value = "{name}: {index}".format(
            name=name,
            index=str(index or '?')
        )

def symbolic_matrix_factor(M):
    import IPython.display
    from IPython.display import display
    import sympy as _sp

    # Factorize each element of the matrix
    factored_M = M.applyfunc(_sp.factor)
    
    # Flatten the matrix into a list
    elements = factored_M.tolist()
    elements = flatten(elements)


    # Remove the zeros from the list
    elements = [elem for elem in elements if elem != 0]

    # Remove duplicates
    elements = list(set(elements))

    # Extract the common factor
    common_factor = 1
    
    while True:
        if len(elements) == 0:
            break
        try:
            current_gcd = _sp.gcd(*elements)
            common_factor *= current_gcd
            elements = [_sp.simplify(_elem / current_gcd) for _elem in elements]
            elements = list(set(elements))
            if current_gcd == 1:
                break
        except _sp.GeneratorsError:
            break

    # Simplify the matrix by dividing by the common factor
    simplified_M = _sp.simplify(factored_M / common_factor)

    # Get LaTex representation
    simplified_M_latex = _sp.latex(simplified_M)
    common_factor_latex = _sp.latex(common_factor)
    
    # check if common factor is a term where the highest level operation is a sum or a difference.
    # if so, we need to wrap it in parentheses
    if common_factor.is_Add:
        common_factor_latex = '\\left(' + common_factor_latex + '\\right)'
    
    latex_expression = '\\begin{equation*}' + common_factor_latex + '\\, \\cdot \\,' + simplified_M_latex + '\\end{equation*}'


    display(IPython.display.Latex(latex_expression))

        
def symbolic_replace_all(expr, value, *keep):
    import sympy as _sp
    subs_dict = {var: value for var in expr.free_symbols if var not in keep}
    expr_mod = expr.subs(subs_dict)
    expr_mod = expr_mod.subs({_sp.physics.quantum.constants.hbar:value})
    return expr_mod