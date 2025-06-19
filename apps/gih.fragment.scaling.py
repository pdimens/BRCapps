import marimo

__generated_with = "0.14.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
    # BRC Fragment Analysis Scaled Concentrations
    This worksheet will scale the concentration of your DNA samples based on the proportion of representation of your target fragment interval as determined by fragment analysis. In other words, given the fragment analysis results, your target interval, and the original concentrations of your samples (or pools), this worksheet will calculate the "effective" concentration of the DNA you want to sequence.
    """
    )
    return


@app.cell(hide_code=True)
def _():
    import io
    import os
    from pathlib import Path
    import marimo as mo
    import pandas as pd

    example_table = mo.md("""
    | Well | Sample ID | Range             | ng/uL  | % Total | nmole/L | Avg. Size | %CV   | Size Threshold (b.p.) | DQN |
    |:-----|:----------|:------------------|:-------|:--------|:--------|:----------|:------|:----------------------|:----|
    | F2   | sample_1  | 10 bp to 100 bp   | 0.4294 | 6.8     | 14.9616 | 47        | 21.23 | 300                   | 8.6 |
    | F2   | sample_1  | 100 bp to 450 bp  | 0.0312 | 0.5     | 0.2444  | 210       | 7.45  | 300                   | 8.6 |
    | F2   | sample_1  | 450 bp to 800 bp  | 2.4305 | 38.5    | 6.9599  | 575       | 16.43 | 300                   | 8.6 |
    | F2   | sample_1  | 800 bp to 5500 bp | 0.6945 | 11.0    | 0.9207  | 1242      | 49.03 | 300                   | 8.6 |
    | F3   | sample_2  | 10 bp to 100 bp   | 0.4318 | 4.8     | 15.2126 | 46        | 24.69 | 300                   | 8.8 |
    | F3   | sample_2  | 100 bp to 450 bp  | 3.579  | 39.4    | 16.2973 | 361       | 14.84 | 300                   | 8.8 |
    | F3   | sample_2  | 450 bp to 800 bp  | 3.5493 | 39.1    | 10.0215 | 583       | 16.55 | 300                   | 8.8 |
    | F3   | sample_2  | 800 bp to 5500 bp | 1.4898 | 16.4    | 1.9335  | 1268      | 44.82 | 300                   | 8.8 |
    """)
    return example_table, io, mo, pd


@app.cell(hide_code=True)
def _(example_table, mo):
    mo.md(
        f"""
    ## Import Data
    Use the file importer below to load the data from the CSV file that was provided to you by the BRC from the fragment analyzer. Then, choose which row (within sample rows) has the fragment size interval you are interested in performing a scaled concentration correction on. In the example table below, if you were interested in 450bp to 800bp, that would be row 3 of `sample_1`. This value is expected to be consistent across all samples, meaning interval 450-800 would also be the 3rd row for `sample_2`, etc.

    {mo.accordion({"Example CSV File": example_table})}
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    file_import = mo.ui.file(
        kind="area",
        filetypes = [".csv", ".CSV"],
        label = "Drag and drop the fragment analysis CSV file here, or click to open file browser"
    )
    rownumber = mo.ui.number(start=1, stop=10, label="Row number of target interval within each sample")
    mo.vstack([file_import, rownumber], heights= [.25, .75], align='center', justify='center')
    return file_import, rownumber


@app.cell(hide_code=True)
def _(file_import, io, mo, pd, rownumber):
    wait_text = """
    /// admonition| Input file required.

    Cannot proceed without a CSV being uploaded
    ///
    """
    mo.stop(not file_import.value, mo.md(wait_text))

    contents = io.BytesIO(file_import.value[0].contents)
    df = pd.read_csv(contents)
    interval = df['Range'][rownumber.value - 1]
    df
    return df, interval


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
    ## Sample Concentrations
    Using the table below, please include the concentrations, in ng/uL (nanograms per microliter), for each sample.
    """
    )
    return


@app.cell(hide_code=True)
def _(df, mo):
    samples = df[['Well', 'Sample ID']].drop_duplicates().reset_index(drop=True)
    samples['concentration (ng/uL)'] = [0.0] * len(samples)
    editor = mo.ui.data_editor(samples)
    editor
    return (editor,)


@app.cell(hide_code=True)
def _(pd):
    def process_sample(group, _rownum):
        target_row = group.iloc[_rownum.value -1]
        corrected_sum = group['ng/uL'].iloc[1:].sum()
        target_conc = target_row['ng/uL']
        target_size = target_row['Avg. Size']
        sample_id = target_row['Sample ID']
        corrected_smear = target_conc/corrected_sum
        quant = target_row['concentration (ng/uL)']
        return pd.Series(
            {
            'Sample ID': sample_id,
            'Avg.Size': target_size,
            'Window ng/uL': round(target_conc,3),
            'Corrected Smear': round(corrected_smear,3),
            'Sample ng/uL': quant,
            'Corrected ng/uL': round(quant * corrected_smear,3)
        })

    return (process_sample,)


@app.cell(hide_code=True)
def _(interval, mo):
    mo.md(
        f"""
    ## Results

    Below is the table of the adjusted concentrations for the fragment interval you are interested in ({interval}). If this is the wrong interval, please change the row number at the top of this page. The numbers below are rounded to 3 decimal places.
    """
    )
    return


@app.cell(hide_code=True)
def _(df, editor, process_sample, rownumber):
    df_with_conc = df.merge(editor.value, on=['Well','Sample ID'], how="left")
    df_with_conc.groupby('Well').apply(lambda group: process_sample(group, rownumber), include_groups=False)
    return


if __name__ == "__main__":
    app.run()
