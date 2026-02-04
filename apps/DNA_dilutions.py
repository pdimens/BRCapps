import marimo

__generated_with = "0.19.7"
app = marimo.App(width="full", app_title="DNA Dilution Calculator")


@app.cell
def _():
    import copy
    import io
    import marimo as mo
    import numpy as np
    import pandas as pd
    import re
    import string
    return copy, io, mo, pd


@app.cell
def _(mo):
    # define UI elements

    file_import = mo.ui.file(
        kind="area",
        filetypes = [".csv", ".CSV", ".Csv", ".txt", ".TXT", ".tsv", ".TSV"]
    )

    start_ul = mo.ui.slider(
        value = 1,
        start = 0.5,
        stop = 200,
        step = 0.5,
        include_input = True,
        full_width = True,
        label = "Microliters (ul) of input DNA"
    )

    target_ng = mo.ui.slider(
        value = 5,
        start = 0.1,
        stop = 10,
        step = 0.1,
        include_input = True,
        full_width = True,
        label = "Target concentration (ng/ul) of DNA"
    )

    headers = mo.ui.switch(value= True, label = "has headers")

    example_table = mo.md("""
    | Well | Sample    | ng/uL  |
    |:-----|:----------|:-------|
    | A1   | sample_1  | 0.4294 |
    | A2   | sample_2  | 0.0312 |
    | A3   | sample_3  | 2.4305 |
    | A4   | sample_4  | 0.6945 |
    """)

    example_file = mo.md("""
    |Well,Sample,ng/uL|
    |:----|
    |A1,sample_1,0.4294|
    |A2,sample_2,0.0312|
    |A3,sample_3,2.4305|
    |A4,sample_4,0.6945|
    """
    )
    return example_file, file_import, headers, start_ul, target_ng


@app.cell
def _(copy, df, pd, start_ul, target_ng):
    # define functions
    def print_upto(x, n) -> list[str]:
        '''Given a stringio bytes object `x`, split out lines and return up to `n` lines'''
        _lines = x.decode("utf-8").split('\n')
        return _lines[:min(n, len(_lines))]

    def parse_well(well):
        '''Function to parse well notation (handles both A1 and A01 as column 1)'''
        row = well[0]  # First character is the row (A-H)
        col = str(int(well[1:]))  # Remaining characters are the column number
        return row, col

    def style_cell(rowId, columnName, value):
        '''Only style cells in the diluent column where value < 0'''
        if columnName == 'diluent (ul)':
            if value < 0:
                return {
                    "backgroundColor": "lightcoral",
                    "color": "darkred",
                    "fontWeight": "bold"
                }
            elif  value > (190 - start_ul.value):
                return {
                    "backgroundColor": "orange",
                    "color": "brown",
                    "fontWeight": "bold"
                }
        return {}

    def style_well(rowId, columnName, value):
        '''Only style cells where value is <0 or >190'''
        if columnName != 'row':
            if value < 0:
                return {
                    "backgroundColor": "lightcoral",
                    "color": "darkred",
                    "fontWeight": "bold"
                }
            elif  value > (190 - start_ul.value):
                return {
                    "backgroundColor": "orange",
                    "color": "brown",
                    "fontWeight": "bold"
                }
        return {}

    def recalc_dilution(updated_df):
        '''callback function for editable dataframe'''
        updated_df['diluent (ul)'] = round(((updated_df['ng/ul'] * updated_df['ul DNA']) / target_ng.value) - updated_df['ul DNA'], 2)

    def table_to_plate(input_table) -> pd.DataFrame:
        '''convert long-form table into 96-well plate format dataframe'''
        df_copy = copy.copy(input_table)
        # Extract row and column from well notation
        df_copy['row'] = df['well'].apply(lambda x: parse_well(x)[0])
        df_copy['col'] = df['well'].apply(lambda x: parse_well(x)[1])

        # Create the transposed table
        plate_layout = df_copy.pivot(index='row', columns='col', values='diluent (ul)')

        # Ensure all rows A-H and columns 1-12 exist (fill missing with 0)
        all_rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        all_cols = [str(i) for i in range(1, 13)]

        return plate_layout.reindex(index=all_rows, columns=all_cols).fillna(0)
    return recalc_dilution, style_cell, style_well, table_to_plate


@app.cell
def _(file_import, headers, mo, start_ul, target_ng):
    mo.sidebar(
        [
            mo.md(
                '# DNA Dilutions Calculator'
                '\nThis worksheet takes your DNA concentrations and formats the water volumes for a dilution plate.'
            ),
            mo.vstack([
              file_import,
              headers,
              start_ul,
              target_ng
            ], align = 'center'),
        ],
        footer = mo.md('<img src="public/gih_logo.png" width="200" />\n\nMade with ❤️ for 🧬')
    )
    return


@app.cell
def _(mo):
    mo.md(f"""
    ## Import Data
    Use the file importer in the left sidebar to load the data from a comma or whitespace delimited (CSV/TSV) file. At minimum, the file needs a column of well ID's (e.g. `A2`), sample names, and concentration (ng/ul). Wells can be in 1 or 2 digit form (e.g. `A1` and `A01` are both valid).
    """)
    return


@app.cell
def _(example_file, file_import, headers, io, mo, pd):
    mo.stop(
        not file_import.value,
        mo.md(
            f"""
    /// admonition| Input file required

    /// details | 🔍︎ Example input
    {example_file}
    ///
    
    ///
            """
        )
    )

    try:
        df = pd.read_table(io.BytesIO(file_import.value[0].contents), header= 0 if headers.value else None, engine='python', sep=None)
        df.columns = [i.lower() for i in df.columns]
    except Exception:
        is_err = True
        err_md = mo.md(f"""
        /// error| Invalid input file

        The `pandas` parser failed to read the input file. Please check that it conforms to a tabular comma or whitespace (tab/space) delimited file. The first 200 bytes of the uploaded file:

        ///
        """)
        err = mo.vstack([err_md, mo.md(f"```\n{file_import.value[0].contents[:200]}\n```")])
        mo.stop(is_err, err)

    colnames = mo.ui.array([mo.ui.dropdown(options=list(df.columns), label=f"{i} column") for i in ['well','sample ID', 'concentration']])

    mo.vstack([
        mo.md("Please select the columns in your input file to map as well, sample ID, and concentration columns"),
        mo.ui.table(df.head(3), selection = None, show_column_summaries=False, show_data_types=False, show_download = False),
        mo.hstack(colnames)
    ])
    return colnames, df


@app.cell
def _(colnames, df, mo, start_ul, target_ng):
    mo.stop(
        any([not i.value for i in colnames]),
        mo.md(f"""
            /// admonition| Map input columns to continue

            ///
            """)
    )

    df_clean = df[[i.value for i in colnames]]
    df_clean.columns = ['well', 'sample', 'ng/ul']
    df_clean['ul DNA'] = start_ul.value
    df_clean['diluent (ul)'] = round(((df_clean['ng/ul'] * df_clean['ul DNA']) / target_ng.value) - df_clean['ul DNA'], 2)
    return (df_clean,)


@app.cell
def _(colnames, df_clean, file_import, mo, recalc_dilution):
    mo.stop(not file_import.value or (any([not i.value for i in colnames])))
    input_text = mo.md(r"""
    ## Input Volume
    Use the slider on the left to set a fixed input volume or the table editor in the dropdown below to modify the input DNA volume  per sample. The changes will be reflected in the final output tables.
    """)

    editor = mo.ui.data_editor(df_clean, editable_columns=['ul DNA'], on_change=recalc_dilution)
    mo.vstack([
        input_text,
        mo.accordion({"Variable Volumes" : editor})
    ])

    return (editor,)


@app.cell
def _(editor, mo, style_cell):
    output_table = mo.ui.table(
        editor.value,
        show_data_types=False,
        pagination = False,
        show_column_summaries=False,
        show_download=True,
        style_cell = style_cell,
        label = "Dilution Values"
    )
    return (output_table,)


@app.cell
def _(colnames, file_import, mo):
    mo.stop(not file_import.value or (any([not i.value for i in colnames])))
    mo.md(r"""
    ## Resulting DNA Dilutions
    Negative values for `diluent (ul)` are highlighted in red and indicate that the sample's original concentration is below your desired diluted concentration. Values for `diluent (ul)` highlighted in orange indicate a final volume greater than **190ul**, which would be the maximum volume a standard 96-well microplate can hold without issues.
    """)
    return


@app.cell
def _(editor, mo, start_ul, style_well, table_to_plate):
    plate_fmt = mo.ui.table(
        table_to_plate(editor.value),
        show_data_types=False,
        selection = None,
        pagination = False,
        show_column_summaries=False,
        show_download = True,
        style_cell = style_well,
        label = f"Microliters (**ul**) of diluent to add to **{start_ul.value}ul** of DNA (unless manually modified)"
    )
    return (plate_fmt,)


@app.cell
def _(mo, pd):
    def download_mantis(input: pd.DataFrame):
        all_cols = [str(i) for i in range(1, 13)]
        res = "," + ",".join(all_cols) + "\n"
        for r in input.itertuples():
            res += ",".join([str(i) for i in r]) + "\n"
        return mo.download(
        data=res.encode("utf-8"),
        filename="manits.dilution.txt",
        mimetype="text/plain",
        label="Download as Mantis Configuration",
    )
    return (download_mantis,)


@app.cell
def _(download_mantis, mo, output_table, plate_fmt):
    mo.ui.tabs(
        {"Plate View": mo.vstack([download_mantis(plate_fmt.data), plate_fmt]), "Table View": output_table},
        label = "Dilution Values", lazy = True
    )
    return


if __name__ == "__main__":
    app.run()
