# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "pandas==2.3.3",
#     "pygnome==0.3.0",
# ]
# ///

import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(r"""
    # GTF File Summarizer
    """)
    return


@app.cell
def _():
    import io
    from itertools import batched
    import gzip
    import marimo as mo
    import os
    import pandas as pd
    return batched, gzip, io, mo, pd


@app.cell
def _(mo):
    file_import = mo.ui.file(
        kind="area",
        filetypes = [".gtf", ".GTF", ".gtf.gz", ".GTF.gz", ".GTF.GZ", ".gff", ".GFF", ".gff.gz", ".GFF.gz", ".GFF.GZ"],
        label = "Drag and drop the GTF file here, or click to open file browser. The maximum file size is 2GB. If your file is larger than 2GB, you can try to gz-compress the file to shrink it.",
        max_size =  2000000000
    )
    file_import
    return (file_import,)


@app.cell
def _(file_import, gzip, io, mo, pd):
    wait_text = """
    /// admonition| Input file required.

    Cannot proceed until a file is uploaded
    ///
    """
    mo.stop(not file_import.value, mo.md(wait_text))

    def process_attributes(x):
        _semicolon_split = x.strip(";").replace("\"",'').split(";")
        d = {}
        for i in _semicolon_split:
            if i.strip():
                #print(i.split())
                keyval = i.strip().split()
                if len(keyval) > 1:
                    d[keyval[0].strip()] = keyval[1].strip()
        return d

    if file_import.name().lower().endswith(".gz"):
        infile = gzip.open(io.BytesIO(file_import.contents()), "r")
    else:
        infile = io.BytesIO(file_import.contents())

    df = pd.read_csv(
        infile,
        delimiter="\t",
        header = None,
        comment = "#",
        skip_blank_lines=True,
        names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"],
        converters= {
            "attribute" : process_attributes
        }
    )

    df.insert(8, 'gene id', [i["gene_id"] for i in df['attribute']])
    for j in df['attribute']:
        del j["gene_id"]

    _generows = (df['feature'] == 'gene').sum()
    mdtext = f"""
    Input file `{file_import.value[0].name}` has **{_generows} gene feature rows**.
    Showing the first 5 rows.
    """
    df.head()
    mo.vstack(
        [mo.md(mdtext), df.head()]
    )
    return (df,)


@app.cell
def _(df, mo):
    allkey = set()
    for a1 in df['attribute']:
        [allkey.add(str(a2)) for a2 in a1]
    allkeys = sorted(allkey)
    switches = mo.ui.array([mo.ui.switch(label=b1) for b1 in allkeys])
    return allkeys, switches


@app.cell
def _(batched, mo, selected_attributes, switches):
    keycols = ", ".join(selected_attributes)
    mo.vstack(
        [
        mo.md(f"## Summary Table\nThe table below summarizes across unique `gene_id` values and reports the `start`/`end` positions that corrspond to the information in the `gene` rows (of the `feature` column). It includes columns consolidating unique values for the attribute names selected by switching the toggles you see below."),
        mo.hstack([mo.vstack(z, gap = 0.05) for z in batched(switches, 10)])
        ]
    )
    return


@app.cell
def _(allkeys, df, mo, pd, switches):
    selected_attributes = [i for idx,i in enumerate(allkeys) if switches.value[idx]]

    generows = df[df['feature'] == 'gene'].iloc[:, [8,0,3,4,6]]

    user_defined_cols = df.groupby('gene id').apply(
        lambda x: pd.Series({
            key: ';'.join(
                sorted(set(
                    str(d[key]) for d in x['attribute'] if key in d
                )
            ))
            for key in selected_attributes
        }),
        include_groups=False
    ).reset_index()

    mo.ui.table(
        generows.merge(user_defined_cols, on='gene id', how='left'),
        show_column_summaries=False
    )
    return (selected_attributes,)


if __name__ == "__main__":
    app.run()
