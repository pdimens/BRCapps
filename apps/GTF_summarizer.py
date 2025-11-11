# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "pandas==2.3.3",
#     "pygnome==0.3.0",
# ]
# ///

import marimo

__generated_with = "0.17.7"
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
    import gzip
    import marimo as mo
    import os
    import pandas as pd
    return gzip, io, mo, pd


@app.cell
def _(mo):
    file_import = mo.ui.file(
        kind="area",
        filetypes = [".gtf", ".GTF", ".gtf.gz", ".GTF.gz", ".GTF.GZ"],
        label = "Drag and drop the GTF file here, or click to open file browser"
    )
    file_import
    return (file_import,)


@app.cell
def _(file_import, gzip, io, mo, pd):
    wait_text = """
    /// admonition| Input file required.

    Cannot proceed until a GTF file is uploaded
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

    if file_import.value[0].name.lower().endswith(".gz"):
        df = pd.read_csv(
            gzip.open(io.BytesIO(file_import.contents()), "r"),
            delimiter="\t",
            header = None,
            comment = "#",
            skip_blank_lines=True,
            names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"],
            converters= {
                "attribute" : process_attributes
            }
        )

    else:
        df = pd.read_csv(
            io.BytesIO(file_import.contents()),
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
    allkeys = set()

    for i in df['attribute']:
        [allkeys.add(str(j)) for j in i]

    multiselect = mo.ui.multiselect(
        options = sorted(allkeys),
        label = "Choose which attribute keys to retain as columns below. Rows with the same `gene_id` that have repeated attribute keys across their rows will have their unique values concatenated using a semicolon `;`",
        full_width=True
    )

    multiselect
    return (multiselect,)


@app.cell
def _(mo, multiselect):
    keycols = ", ".join(multiselect.value)
    mo.md(f"# Summary Table\nThis table summarizes across unique `gene_id` values and reports the `start`/`end` positions that corrspond to the information in the `gene` rows (of the `feature` column). It includes columns consolidating unique values for the values in the `attribute` keys:\n**{keycols}**")
    return


@app.cell
def _(df, mo, multiselect, pd):
    generows = df[df['feature'] == 'gene'].iloc[:, [8,0,3,4,6]]

    user_defined_cols = df.groupby('gene id').apply(
        lambda x: pd.Series({
            key: ';'.join(
                sorted(set(
                    str(d[key]) for d in x['attribute'] if key in d
                )
            ))
            for key in multiselect.value
        }),
        include_groups=False
    ).reset_index()

    mo.ui.table(
        generows.merge(user_defined_cols, on='gene id', how='left'),
        show_column_summaries=False
    )
    return


if __name__ == "__main__":
    app.run()
