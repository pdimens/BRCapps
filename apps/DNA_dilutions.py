import marimo

__generated_with = "0.19.7"
app = marimo.App(width="full", app_title="DNA Dilution Calculator")


@app.cell
def _():
    import copy
    import io
    import marimo as mo
    from math import modf
    import numpy as np
    import pandas as pd
    import re
    import string
    return copy, io, mo, modf, np, pd


@app.cell
def _(mo):
    # define UI elements

    file_import = mo.ui.file(
        kind="area",
        filetypes = [".csv", ".CSV", ".Csv", ".txt", ".TXT", ".tsv", ".TSV"]
    )

    start_ul = mo.ui.slider(
        value = 1,
        start = 0.1,
        stop = 25,
        step = 0.1,
        include_input = True,
        full_width = True,
        label = "Input volume (ul) of stock DNA"
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

    target_ul = mo.ui.slider(
        value = 5,
        start = 1,
        stop = 20,
        step = 0.1,
        include_input = True,
        full_width = True,
        label = "Min. required ul of diluted DNA"
    )


    headers = mo.ui.switch(value= True, label = "file has column headers")

    example_file = mo.md("""
    |Well,Sample,ng/uL|
    |:----|
    |A1,sample_1,0.4294|
    |A2,sample_2,0.0312|
    |A3,sample_3,2.4305|
    |A4,sample_4,0.6945|
    """
    )
    return example_file, file_import, headers, start_ul, target_ng, target_ul


@app.cell
def _(copy, df, modf, pd, start_ul, target_ng, target_ul):
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
                    "backgroundColor": "black",
                    "color": "white",
                    "fontWeight": "bold"
                }
            elif (value + start_ul.value) < target_ul.value:
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
                    "backgroundColor": "black",
                    "color": "white",
                    "fontWeight": "bold"
                }
            elif (value + start_ul.value) < target_ul.value:
                return {
                    "backgroundColor": "orange",
                    "color": "brown",
                    "fontWeight": "bold"
                }
        return {}

    def recalc_dilution(updated_df):
        '''callback function for editable dataframe'''
        updated_df['diluent (ul)'] = round(((updated_df['ng/ul'] * updated_df['ul DNA']) / target_ng.value) - updated_df['ul DNA'], 2)
        updated_df['total volume (ul)'] = round(updated_df['diluent (ul)'] + start_ul.value, 1)

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

    def mantis_round(x):
        '''
        Given a input integer or float, return a list of `[integer, decimal]`, e.g. `4.2` -> `[4,0.2]`.
        Returns `[0,0]` if `x` < 0
        '''
        if x < 0 or x > (190 - start_ul.value):
            return [0,0]
        b,a = modf(x)
        _b = round(b, 1) if b > 0 else 0
        return [int(a), _b]
    return (
        mantis_round,
        recalc_dilution,
        style_cell,
        style_well,
        table_to_plate,
    )


@app.cell
def _(file_import, headers, mo, start_ul, target_ng, target_ul):
    mo.sidebar(
        [
            mo.md(
                '# DNA Dilutions Calculator'
                '\nThis worksheet takes your DNA concentrations and formats the water volumes for a dilution plate.'
            ),
            file_import,
            headers,
            start_ul,
            target_ng,
            target_ul
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
    df_clean['diluent (ul)'] = round(((df_clean['ng/ul'] * df_clean['ul DNA']) / target_ng.value) - df_clean['ul DNA'], 1)
    df_clean['total volume (ul)'] = round(df_clean['diluent (ul)'] + start_ul.value, 1)
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
def _(colnames, file_import, mo, target_ng, target_ul):
    mo.stop(not file_import.value or (any([not i.value for i in colnames])))
    mo.md(fr"""
    ## Resulting DNA Dilutions
    - <strong style="background-color:lightcoral;color:darkred;">Red</strong> highlighted wells have a starting concentration below the target diluted concentration of **{target_ng.value}ul**
    - <strong style="background-color:orange;color:brown;">Orange</strong> highlighted wells have a final volume below the target minimum volume of **{target_ul.value}ul**
    - <strong style="background-color:black;color:white;">Black</strong> highlighted wells have a final volume greater than **190ul**, which would be the maximum volume a standard 96-well microplate can hold without issues
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
def _(mantis_round, mo, np, pd):
    def format_mantis(arr):
        res = []
        for row in np.array_split(arr, 8):
            res.append("\t".join(str(i) for i in row))
        return "\n".join(res)

    def download_mantis(input: pd.DataFrame):
        outfile = ["[ Version: 5 ]\nAB0800-96well-on-silver-metal-base.pd.txt"]
        high_reagent = "Wash + T HV		Normal"
        low_reagent = "Wash + T LV		Normal"
        high_vec = []
        high_tracks = []
        res_low = ""
        for r in input.itertuples(index = False):
            low_chip = []
            high_chip = []
            for i in r:
                vals = mantis_round(i)
                high_chip.append(vals[0])
                low_chip.append(str(vals[1]))
            high_vec += high_chip  
            res_low  += "\t".join(low_chip) + "\n"
        track_vol = 0
        cutoff = 850
        allwells = [0] * 96
        for i,well_vol in enumerate(high_vec):
            track_vol += well_vol
            if track_vol >= cutoff or i == 95:
                high_tracks.append(format_mantis(allwells))
                # reset
                track_vol = well_vol
                allwells = [0] * 96
            allwells[i] = well_vol
        outfile.append("5\t0\t\t" + "\t\t".join(["U"] * len(high_tracks)) + "\t")
        outfile.append("1")
        outfile.append("5\t0\t\t" + "\t\t".join(["0"] * len(high_tracks)) + "\t")
        for i in high_tracks:
            outfile.append(high_reagent)
            outfile.append("Well\t1")
            outfile.append(i)
        outfile.append(low_reagent)
        outfile.append("Well\t1")
        outfile.append(res_low)
        return mo.download(
            data=("\n".join(outfile)).encode("utf-8"),
            filename="manits.dilution.dl.txt",
            mimetype="text/plain",
            label="Download Mantis config"
    )
    return (download_mantis,)


@app.cell
def _(download_mantis, mo, output_table, plate_fmt):
    mantis_dl = mo.vstack([
        download_mantis(plate_fmt.data),
        mo.md("Use the button above to download the dilution plate as a Mantis configuration file. The configuration sets up a series of steps for the high-volume chip to dispense the diluent up to the integer value of the volume (e.g. `40` of `40.1` ul). The high volume step will be separated into smaller steps to make sure the total dispensed volume does not exceed 850ul so that the diluent can be topped off and proceed to the next series of wells. The final step uses the low-volume chip to dispense the diluent for the remaining <1ul volume (e.g. `0.1` of `40.1` ul). Negative values and values with total expected volumes > `190` will be replaced with `0`.")
    ])

    mo.ui.tabs(
        {"Plate View": mo.vstack([plate_fmt, mantis_dl]), 
         "Table View": output_table
        },
        label = "Dilution Values", lazy = True
    )
    return


if __name__ == "__main__":
    app.run()
