# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "pandas==2.3.3",
# ]
# ///

import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import gzip
    import io
    from itertools import batched
    import marimo as mo
    import os
    import pandas as pd
    import re
    pd.options.mode.copy_on_write = True
    return batched, gzip, io, mo, pd, re


@app.cell
def _(mo):
    file_import = mo.ui.file(
        kind="area",
        filetypes = [".gtf", ".GTF", ".gtf.gz", ".GTF.gz", ".GTF.GZ", ".gff", ".GFF", ".gff.gz", ".GFF.gz", ".GFF.GZ"],
        label = "Drag and drop the GTF file here, or click to open file browser.",
        max_size =  2000000000
    )
    return (file_import,)


@app.cell
def _(file_import, mo):
    mo.sidebar([mo.md('# GTF Summarizer\nThe maximum file size is 2GB. If your file is larger than 2GB, you can try to gz-compress the file to shrink it.'), file_import], footer = mo.md('<img src="public/gih_logo.png" width="200" />\n\nMade with ❤️ for 🧬'))
    return


@app.cell
def _(re):
    def parse_attributes(text):
        """Process the attributes column of a GTF file and correctly assign key-value pairs"""
        result = {}
        # Find all key followed by one or more quoted values
        pattern = r'(\w+)\s+((?:"[^"]+"\s*)+)'

        for match in re.finditer(pattern, text):
            key = match.group(1)
            values_str = match.group(2)

            # Extract all quoted values for this key
            values = re.findall(r'"([^"]+)"', values_str)

            # Store as single value or list
            result[key] = values[0] if len(values) == 1 else ", ".join(values)

        return result
    return (parse_attributes,)


@app.cell
def _(file_import, gzip, io, mo, parse_attributes, pd):
    wait_text = """/// admonition| Input file required.\n\nCannot proceed until a file is uploaded\n///"""
    mo.stop(not file_import.value, mo.md(wait_text))

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
            "attribute" : parse_attributes
        }
    )

    df.insert(8, 'gene id', [i["gene_id"] for i in df['attribute']])
    for j in df['attribute']:
        del j["gene_id"]

    _generows = (df['feature'] == 'gene').sum()
    mdtext = f"Input file `{file_import.value[0].name}` has **{_generows} gene feature rows**. Showing the first 5 rows."""
    mo.ui.table(df.head(), show_data_types=False, label = mdtext)
    return (df,)


@app.cell
def _(df, mo):
    allkey = set()
    for a1 in df['attribute']:
        [allkey.add(str(a2)) for a2 in a1]
    allkeys = sorted(allkey)
    switches = mo.ui.array([mo.ui.switch(label=b1) for b1 in allkeys])
    button = mo.ui.run_button(label = "Generate output table")
    return allkeys, button, switches


@app.cell
def _(batched, mo, switches):
    mo.vstack(
        [
        mo.md(f"## Summary Table\nThe table below summarizes across unique `gene_id` values and reports the `start`/`end` positions that corrspond to the information in the `gene` rows (of the `feature` column). It includes columns consolidating unique values for the attribute names selected by switching the toggles you see below. Once you select the attributes you want, press the \"**Generate output table**\" button below to create the table (this is done to avoid recomputing the table each time an attribute is selected)."),
        mo.hstack([mo.vstack(z, gap = 0.05) for z in batched(switches, len(switches) // 3 + 1)])
        ]
    )
    return


@app.cell
def _(allkeys, button, df, mo, pd, switches):
    mo.stop(not button.value, output = button.center())
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
        label = "Summary Table",
        page_size = 20,
        show_column_summaries=False,
        show_data_types=False
    )
    return


if __name__ == "__main__":
    app.run()
