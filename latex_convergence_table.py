#==============================================================================
# Available norms. Some of these may be missing, but they always follow this
# order when printed
available_norms = ['l2', 'h1', 'h2']

# How to print the various norms
norm_labels = {'l2': '$L^2$',
               'h1': '$H^1$',
               'h2': '$H^2$'}

#==============================================================================
def print_ncell(n):
    return '${:d}$'.format(n)

def print_error(e):
    float_str = "{:8.2e}".format(e)
    mantissa, exponent = float_str.split('e')
    return r'${0} \times 10^{{{1}}}$'.format(mantissa, int(exponent))

def print_order(o):
    if o is None:
        return ' -- '
    else:
        return '${:3.2f}$'.format(o)

#==============================================================================
def latex_tabulate(table, headers, alignment):

    lines = []

    line = r'\begin{tabular}{%s}' % alignment
    lines.append(line)

    line = r'\hline'
    lines.append(line)

    for header in headers:
        line = ' & '.join( header ) + r' \\'
        lines.append(line)

    line = r'\hline'
    lines.append(line)

    for entries in table:
        line = ' & '.join( entries ) + r' \\'
        lines.append( line )

    line = r'\hline'
    lines.append(line)

    line = r'\end{tabular}'
    lines.append(line)

    return '\n'.join(lines)

#==============================================================================
def convergence_table(n_list, **errors):

    import math
    # Check that all norms are recognized
    for norm in errors.keys():
        if norm not in available_norms:
            raise Warning('Cannot recognize norm', norm)

    # Store norm names and LaTeX labels into lists
    used_norms  = []
    used_labels = []
    for norm in available_norms:
        if norm in errors.keys():
            used_norms .append(norm)
            used_labels.append(norm_labels[norm])

    # Compute convergence order between two data points
    order = lambda n1, n2, e1, e2: math.log(e2/e1) / math.log(n1/n2)

    # Compute dictionary with convergence orders for all norms
    orders = {}
    for norm, err_list in errors.items():
        orders[norm] = [None] + [order(n1, n2, e1, e2) for n1, n2, e1, e2 in
                zip(n_list[1:], n_list[:-1], err_list[1:], err_list[:-1])]

    # Data organized by columns
    columns = [n_list] + [d[norm] for norm in used_norms for d in (errors, orders)]

    # Format line entries
    def print_line(n, *error_order):

        errors = error_order[0::2]
        orders = error_order[1::2]

        line = [print_ncell(n)]
        for e, o in zip(errors, orders):
            line += [print_error(e), print_order(o)]

        return line

    # Create table:
    #  . Table is list of lines
    #  . Line is list of entries
    table = [print_line(*entries) for entries in zip(*columns)]

    # Table header
    headers = ['N'] + [r'{0:s} {1:s}'.format(label, what)
                       for label in used_labels
                       for what  in ('error', 'order')]

    # Column alignment
    alignment = 'r' + ' c c' * len(used_norms)

    # LaTeX 'tabulate' object
    txt = latex_tabulate(table, (headers,), alignment)

    return txt

#==============================================================================
def minipage_grid( blocks, subcaptions, ncols ):

    import math

    assert ncols >= 1
    assert len(subcaptions) == len(blocks)

    nblocks = len(blocks)
    nrows   = math.ceil( nblocks / ncols )
    width   = 1.0 / ncols - 0.01
    chunks  = []

    k = -1
    for i in range(nrows):

        # New row of minipages
        for j in range(ncols):

            # Update counter and stop cycle if necessary
            k += 1
            if k == nblocks:
                break

            # Create new minipage
            minipage = [
              r'\begin{minipage}[t]' + r'{{{:.2f}\textwidth}}'.format(width),
              r'\small',
              r'\centering',
              r'\subcaption{{{:s}}}'.format(subcaptions[k]),
              blocks[k],
              r'\end{minipage}'
            ]

            # Store text chunks
            chunks.extend( minipage )

        # Extra spacing between rows
        if i != nrows-1:
            chunks.append(r'\\[1em]')

    # Merge text chunks with 'newline' characters
    txt = '\n'.join(chunks) + '\n'

    return txt

#==============================================================================
def main(*, basename, ncols, varname):

    import re
    import numpy as np
    from pathlib import Path

    ROOT_DIR   = Path(__file__).resolve().parents[0]
    DATA_DIR   = ROOT_DIR
    TABLES_DIR = ROOT_DIR

    pattern = 'errors_' + basename + '_p=([0-9]+).npy'
    regex   = re.compile(pattern)

    paths = [f for f in DATA_DIR.iterdir() if regex.match(f.name)]
    paths.sort()

    if varname:
        key_template = '_'.join(['{norm}', varname])
        tex_filename = '_'.join(['table', basename, varname]) + '.tex'
    else:
        key_template = '_'.join(['{norm}'])
        tex_filename = '_'.join(['table', basename]) + '.tex'

    tables = []
    subcaptions = []

    for path in paths:
        f = np.load(path, allow_pickle=True).item(0)
        kwargs = {}
        keys   = f.keys()
        for norm in available_norms:
            key = key_template.format(norm=norm)
            print(key)
            if key in f:
                kwargs[norm] = f[key]

        table = convergence_table(f['n_list'], **kwargs)

        p = regex.match(path.name).group(1)
        subcaption = 'Degree $p = {}$'.format(p)

        tables     .append( table      )
        subcaptions.append( subcaption )

    # Create LaTeX table and write it to file
    txt  = minipage_grid(tables, subcaptions, ncols)
    path = TABLES_DIR / tex_filename
    with open(path, 'w') as f:
        print(txt, file=f)

#==============================================================================
def parse_input_arguments():

    import argparse

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description     = "Create LaTeX table environment containing " +
                          "a grid of convergence subtables."
    )

    parser.add_argument( 'basename',
        type    = str,
        metavar = 'BASENAME',
        help    = 'Base name of error files'
    )

    parser.add_argument( '-c',
        type    = int,
        default = 2,
        dest    = 'ncols',
        help    = 'Number of columns in grid of subtables'
    )

    parser.add_argument( '-n',
        type    = str,
        default = None,
        dest    = 'varname',
        help    = 'Name of variable of interest (if any)'
    )

    return parser.parse_args()

#==============================================================================
if __name__ == '__main__':
    args = parse_input_arguments()
    main(**vars(args))
