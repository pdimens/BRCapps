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
    import pyarrow
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
    Use the file importer below to load the data from the CSV file that was provided to you by the BRC from the fragment analyzer. Then, choose which row (within sample rows) has the fragment size interval you are interested in performing a scaled concentration correction on. In the example table below, if you were interested in the `Range` 450bp to 800bp, that would be row 3 of `sample_1`. This value is expected to be consistent across all samples, meaning interval 450-800 would also be the 3rd row for `sample_2`, etc.

    {mo.accordion({"Example CSV File": example_table})}
    """
    )
    return


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
            'Corrected ng/uL': round(quant * corrected_smear,3),
            'Est nM': round(quant * corrected_smear,3) / (660 * target_size) * 1000000
        }
    )

    return (process_sample,)


@app.cell(hide_code=True)
def _(mo):
    file_import = mo.ui.file(
        kind="area",
        filetypes = [".csv", ".CSV"],
        label = "Drag and drop the fragment analysis CSV file here, or click to open file browser"
    )
    file_import
    return (file_import,)


@app.cell(hide_code=True)
def _(file_import, io, mo, pd):
    wait_text = """
    /// admonition| Input file required.

    Cannot proceed without a CSV being uploaded
    ///
    """
    mo.stop(not file_import.value, mo.md(wait_text))

    contents = io.BytesIO(file_import.value[0].contents)
    df = pd.read_csv(contents)
    mo.stop(
        sorted(list(df.columns)) != sorted(['Well', 'Sample ID', 'Range', 'ng/uL', '% Total', 'nmole/L','Avg. Size', '%CV', 'Size Threshold (b.p.)', 'DQN']),
        output= mo.md("""
        /// warning| Unrecognized input file

        The input file for this worksheet is expected to have a specific format. It is expected to have these columns, regardless of order:
    
        |Well | Sample ID | Range | ng/uL | % Total | nmole/L | Avg. Size | %CV | Size Threshold (b.p.) | DQN |
        |:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
        ///
        """)
    )
    return (df,)


@app.cell(hide_code=True)
def _(df, mo):
    sample_id = list(set(df['Sample ID']))
    intervals = list(set(df['Range']))
    rownumber = mo.ui.number(start=1, stop=len(intervals), label="Row number of target Range within each sample")
    return (rownumber,)


@app.cell(hide_code=True)
def _(df, mo, rownumber):
    interval = df['Range'][rownumber.value - 1]
    mo.hstack([rownumber, mo.md(f"Range: {interval}")], justify = 'start')
    return (interval,)


@app.cell
def _(df):
    df
    return


@app.cell(hide_code=True)
def _(mo):
    input_header = mo.md(
        """
    ## Sample Concentrations
    Using this interactive form, include the concentrations for each sample in **ng/uL** (nanograms per microliter).
    """
    )
    return (input_header,)


@app.cell(hide_code=True)
def _(df, input_header, mo):
    samples_id = set(df['Sample ID'])
    quants = mo.ui.dictionary(
        dict(zip(samples_id, (mo.ui.text(value = "0", full_width=True,debounce=True) for i in samples_id)))
    )
    mo.vstack([input_header,quants])
    return (quants,)


@app.cell(hide_code=True)
def _(mo, quants):
    err = {}
    for i,j in quants.value.items():
        try:
            float(j)
        except ValueError:
            err[i] = j

    mo.stop(err,
        mo.md(f"""
    /// error| Concentrations must be numeric

    These samples were provided incorrect concentrations:

    {"\n".join(f"{i}: **{j}**" for i,j in err.items())}
    ///""")
    )
    return (err,)


@app.cell(hide_code=True)
def _(err, mo, pd, quants):
    mo.stop(err)
    quants_df = pd.DataFrame(
        {"Sample ID": list(quants.value.keys()), "concentration (ng/uL)": [float(i) for i in quants.value.values()]}
    )
    return (quants_df,)


@app.cell(hide_code=True)
def _(df, err, interval, mo, process_sample, quants_df, rownumber):
    mo.stop(err)
    results_header = mo.md(f"""
    ## Scaled Concentrations

    This table scales the concentrations you input above with the proportion of the sample with the target Range ({interval})
    """)

    df_with_conc = df.merge(
        quants_df,
        on='Sample ID',
        how="left"
    )
    mo.vstack([
        results_header,
        df_with_conc.groupby('Well').apply(lambda group: process_sample(group, rownumber), include_groups=False)
    ])
    return


if __name__ == "__main__":
    app.run()
