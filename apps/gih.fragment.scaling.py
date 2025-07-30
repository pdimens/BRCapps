import marimo

__generated_with = "0.14.13"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
    # BRC Fragment Analysis Scaled Concentrations
    This worksheet will scale the concentration of your DNA samples based on the proportion of representation of your target fragment interval as determined by fragment analysis. In other words, given the fragment analysis results, your target interval, and the original concentrations of your samples (or pools), this worksheet will calculate the "effective" concentration of the DNA you want to sequence.
    """
    )
    return


@app.cell
def _():
    import io
    #import os
    #import pyarrow
    #from pathlib import Path
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


@app.cell
def _(example_table, mo):
    mo.md(
        f"""
    ## Import Data
    Use the file importer below to load the data from the CSV file that was provided to you by the BRC from the fragment analyzer. Then, choose which row (within sample rows) has the fragment size interval you are interested in performing a scaled concentration correction on. In the example table below, if you were interested in the `Range` 450bp to 800bp, that would be row 3 of `sample_1`. This value is expected to be consistent across all samples, meaning interval 450-800 would also be the 3rd row for `sample_2`, etc.

    {mo.accordion({"Example CSV File": example_table})}
    """
    )
    return


@app.cell
def _(pd):
    def process_sample(group, _rownum, picomoles):
        target_row = group.iloc[_rownum.value -1]
        corrected_sum = group['ng/µL'].iloc[1:].sum()
        target_conc = target_row['ng/µL']
        target_size = target_row['Avg. Size']
        sample_id = target_row['Sample ID']
        corrected_smear = target_conc/corrected_sum
        quant = target_row['concentration (ng/µL)']
        nM = (quant * corrected_smear) / (target_size * 660) * 1000000
        if nM > 0:
            vol_to_pool = picomoles / nM
        else:
            vol_to_pool = 0
        ng_primary_lib = quant * vol_to_pool
        return pd.Series(
            {
            'Sample ID': sample_id,
            'Avg.Size': target_size,
            '% of Total Conc.': round(corrected_smear * 100,1),
            'Window ng/µL': round(target_conc,3),
            'Sample ng/µL': quant,
            'Est. nM': round(nM, 2),
            'Corrected ng/µL': round(quant * corrected_smear,3),
            'Volume to Pool': round(vol_to_pool, 2),
            'ng Primary Library': round(ng_primary_lib, 1)
        })

    return (process_sample,)


@app.cell
def _(mo):
    file_import = mo.ui.file(
        kind="area",
        filetypes = [".csv", ".CSV"],
        label = "Drag and drop the fragment analysis CSV file here, or click to open file browser"
    )
    file_import
    return (file_import,)


@app.cell
def _(file_import, io, mo, pd):
    wait_text = """
    /// admonition| Input file required.

    Cannot proceed without a CSV being uploaded
    ///
    """
    mo.stop(not file_import.value, mo.md(wait_text))

    contents = io.BytesIO(file_import.value[0].contents)
    df = pd.read_csv(contents)
    df.rename(columns = {'ng/uL' : 'ng/µL'}, inplace = True)

    mo.stop(
        sorted(list(df.columns)) != sorted(['Well', 'Sample ID', 'Range', 'ng/µL', '% Total', 'nmole/L','Avg. Size', '%CV', 'Size Threshold (b.p.)', 'DQN']),
        output= mo.md("""
        /// warning| Unrecognized input file

        The input file for this worksheet is expected to have a specific format. It is expected to have these columns, regardless of order:

        |Well | Sample ID | Range | ng/µL | % Total | nmole/L | Avg. Size | %CV | Size Threshold (b.p.) | DQN |
        |:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
        ///
        """)
    )
    return (df,)


@app.cell
def _(df):
    sample_id = list(set(df['Sample ID']))
    intervals = list(set(df['Range']))
    return (intervals,)


@app.cell
def _(df):
    df
    return


@app.cell
def _(mo):
    input_header = mo.md("## Sample Concentrations\n\nUsing this interactive form, include the concentrations for each sample in **ng/µL** (nanograms per microliter).")
    return (input_header,)


@app.cell
def _(df, input_header, mo, pd):
    def unique(sequence):
        seen = set()
        return [x for x in sequence if not (x in seen or seen.add(x))]
    samples_id = unique(df['Sample ID'])

    quants_df = mo.ui.data_editor(
        pd.DataFrame({
            'Sample ID' : list(samples_id),
            'concentration (ng/µL)': [0.0 for i in range(len(samples_id))]
        }),
        label = "Add sample concentrations here"
    )

    mo.vstack([input_header,quants_df])
    return (quants_df,)


@app.cell
def _(intervals, mo):
    rownumber = mo.ui.number(start=1, stop=len(intervals), label="Row number of target Range within each sample")
    target_pmol = mo.ui.slider(
        value = 15.0,
        start = 0.1,
        stop = 50.0,
        step = 0.1,
        include_input = True,
        label = "Target **picomoles** for final Libraries"
    )
    return rownumber, target_pmol


@app.cell
def _(df, mo, rownumber, target_pmol):
    interval = df['Range'][rownumber.value - 1]

    mo.vstack([
        mo.hstack([rownumber, mo.md(f"Range: {interval}")], justify = 'start'),
        target_pmol
    ])
    return (interval,)


@app.cell
def _(df, interval, mo, process_sample, quants_df, rownumber, target_pmol):
    df_with_conc = df.merge(
        quants_df.value,
        left_on='Sample ID',
        right_on='Sample ID',
        how="left"
    )

    calc_table = df_with_conc.groupby('Well').apply(lambda group: process_sample(group, rownumber,target_pmol.value), include_groups=False)
    mo.ui.table(
        calc_table,
        pagination = False,
        selection = None,
        show_column_summaries = False,
        show_data_types = False,
        freeze_columns_left = ["Sample ID"],
        label = f"## Scaled Concentrations\n\nThis table scales the concentrations you input above with the proportion of the sample with the target Range ({interval})"
    )
    return (calc_table,)


@app.cell
def _(file_import, mo):
    mo.stop(not file_import.value)
    elution_vol = mo.ui.slider(
        value = 20,
        start = 1,
        stop = 100,
        step = 1,
        include_input = True,
        label = "Volume to elute in (µL)"
    )
    elution_vol
    return (elution_vol,)


@app.cell
def _(calc_table, elution_vol, mo, target_pmol):
    total_vol = calc_table['Volume to Pool'].sum()
    total_ng = calc_table['ng Primary Library'].sum()
    total_pM = target_pmol.value * len(calc_table.index)

    pool_ngul = 0 if total_vol == 0 else round(total_ng / total_vol,1)
    pool_uM = 0 if total_vol == 0 else total_pM / total_vol
    recovery = round(pool_uM * total_vol / elution_vol.value, 2)

    mo.md(f"""
    ## Final Pooling Metrics

    Total Volume of Pool (µL): **{round(total_vol,2)}**

    Total ng Across Pools: **{total_ng}**

    Pool ng/µL : **{pool_ngul}**

    Total pM across intervals: **{total_pM}**

    µM Per Pool across intervals : **{round(pool_uM,1)}**

    µM Per Pool assuming 100% recovery: **{recovery}**
    """)

    return


if __name__ == "__main__":
    app.run()
