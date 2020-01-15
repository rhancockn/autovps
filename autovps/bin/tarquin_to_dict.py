#!/usr/bin/env python
# Converts Tarquin output CSV files to a JSON structure
import pandas as pd
import re


def _df_to_dict(t, prefix='', d=None):
    """Convert a table to a dict.

    Arguments:
        t (pd.DataFrame): table
        prefix (str): string to prefix to each column name
        d (dict): a dictionary to insert keys into
    """

    if d is None:
        d = {}

    if prefix is not '':
        prefix = prefix + '_'

    re_paren = re.compile(r' *\(.*\)')
    for c in t.columns:
        c_name = prefix + re.sub(re_paren, '', c).replace(' ', '_')
        try:
            v = float(t[c][0])
        except ValueError:
            v = t[c][0]
        d[c_name] = v

    return(d)


def tarquin_to_dict(fname):
    """ Return a Tarquin csv file as a dictionary.

    Arguments:
        fname (str): path to Tarquin csv output

    Returns:
        dict (dict)
    """

    # figure out how many columns there are
    n_fields = 0
    with open(fname, 'r') as fp:
        for line in fp:
            n_fields = max(n_fields, len(line.split(',')))

    df = pd.read_csv(fname, header=None, names=range(n_fields),
                     engine='python')
    table_names = ["Signal amplitudes", "CRLBs (standard deviation)",
                   "Fit diagnostics", "Basis shifts and dampings",
                   "Dynamic frequency corrections", "Geometry parameters",
                   "Data parameters"]
    groups = df[0].isin(table_names).cumsum()
    tables = {g.iloc[0, 0]: g.iloc[1:] for k, g in df.groupby(groups)}
    for t_name in tables:
        tables[t_name].columns = tables[t_name].iloc[0]
        tables[t_name] = tables[t_name].reset_index(drop=True)
        tables[t_name] = tables[t_name].reindex(tables[t_name].index.drop(0))
        tables[t_name] = tables[t_name].reset_index(drop=True)
        tables[t_name] = tables[t_name].dropna(how='all', axis=1)

    t_dict = {}
    for t_name in ["Signal amplitudes", "Fit diagnostics"]:
        t_dict = _df_to_dict(tables[t_name], d=t_dict)

    t = tables["CRLBs (standard deviation)"]
    t_dict = _df_to_dict(t, prefix='CRLB', d=t_dict)

    # convert errors to percentages
    crlb_dict = _df_to_dict(t, prefix='CRLB', d={})
    for crlb_metab, crlb_au in crlb_dict.items():
        if 'CRLB_' in crlb_metab:
            metab = '_'.join(crlb_metab.split('_')[1::])
            try:
                t_dict['CRLBPCT_' + metab] = crlb_au / t_dict[metab] * 100.0
            except ZeroDivisionError as e:
                pass

    if "Dynamic frequency corrections" in tables.keys():
        t = tables["Dynamic frequency corrections"]
        delta = pd.to_numeric(t['shift (Hz)'])
        t_dict['fdrift_mean'] = delta.mean()
        t_dict['fdrift_sd'] = delta.std()
        delta_z = (delta - t_dict['fdrift_mean'])/t_dict['fdrift_sd']
        t_dict['fdrift_noutlier'] = sum(delta_z.abs() > 3)
    return(t_dict)

# Run the module as a command line interface
# CSV to JSON

def main():
    import json
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert a Tarquin CSV to a JSON')

    parser.add_argument('--csv', type=str, required=True,
                        help="Path to Tarquin CSV")
    parser.add_argument('--json', type=str, required=True,
                        help="Path to save JSON")
    args = parser.parse_args()

    t_dict = tarquin_to_dict(args.csv)
    fp = open(args.json, 'w+')
    json.dump(t_dict, fp)
    fp.close()


if __name__ == '__main__':
    main()
