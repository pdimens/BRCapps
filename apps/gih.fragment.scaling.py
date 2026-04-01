import marimo

__generated_with = "0.21.1"
app = marimo.App(width="medium", app_title="Fragment Analysis Library Pooling")


@app.cell
def _():
    import io
    import marimo as mo
    import matplotlib as plt
    import numpy as np
    import pandas as pd
    import re
    import statsmodels.formula.api as smf

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
    return example_table, io, mo, np, pd, plt, re


@app.cell
def _(mo):
    file_import = mo.ui.file(
        kind="area",
        filetypes = [".csv", ".CSV", ".Csv"],
        label = "Drag and drop the fragment analysis CSV file here, or click to open file browser"
    )
    return (file_import,)


@app.cell
def _(file_import, mo):
    mo.sidebar(
        [
            mo.md('# Fragment Analysis Scaled Concentrations\nThis worksheet scales the concentration of your samples based on the proportion of representation of your target fragment interval as determined by fragment analysis.'),
            file_import
        ],
        footer = mo.md('<img src="public/gih_logo.png" width="200" />\n\nMade with ❤️ for 🧬')
    )
    return


@app.cell
def _(example_table, mo):
    mo.md(f"""
    ## Import Data
    Use the file importer in the left sidebar to load the data from the CSV file that was provided to you by the BRC from the fragment analyzer. Then, choose which fragment size interval you are interested in performing a scaled concentration correction on.

    {mo.accordion({
        "🔍︎ Example CSV File": example_table,
        "⚠️ Notes and Considerations" : mo.md("**First row skipped**: Be aware that the first interval (sorted alphanumerically) of each sample, usually 10bp-100bp, is skipped in the calculations below.\n\n**Singletons skipped**: Rows with a `Sample ID` that only appears once are removed (e.g. a ladder, samples from a different run).\n\n**Consistency**: The target range is expected to be consistent across all samples.")
    })}
    """)
    return


@app.cell
def _(pd, re):
    def process_sample(group, target_identifier, all_intervals, picomoles):
        target_row = group[group['Range'] == target_identifier].iloc[0]
        corrected_sum = group[group['Range'] != all_intervals[0]]['ng/µL'].sum()
        #group[group['Range'] != '10 to 100 bp']['ng/µL'].sum()
        #corrected_sum = group['ng/µL'].iloc[1:].sum()
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

    def natural_sort(l): 
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(l, key=alphanum_key)

    return natural_sort, process_sample


@app.cell
def _(file_import, io, mo, pd):
    mo.stop(
        not file_import.value,
        mo.md("/// admonition| Input file required\n\nCannot proceed until a CSV file is uploaded\n///")
    )

    df = pd.read_csv(io.BytesIO(file_import.value[0].contents))
    df.rename(columns = {'ng/uL' : 'ng/µL'}, inplace = True)
    df = df[df['Sample ID'].duplicated(keep=False)]
    mo.stop(
        sorted(list(df.columns)) != sorted(['Well', 'Sample ID', 'Range', 'ng/µL', '% Total', 'nmole/L','Avg. Size', '%CV', 'Size Threshold (b.p.)', 'DQN']),
        output= mo.md("""
        /// error| Unrecognized input file

        The input file for this worksheet is expected to have a specific format. It is expected to have these columns, regardless of order:

        |Well | Sample ID | Range | ng/µL | % Total | nmole/L | Avg. Size | %CV | Size Threshold (b.p.) | DQN |
        |:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
        ///
        """)
    )
    return (df,)


@app.cell
def _(df, natural_sort):
    sample_id = list(set(df['Sample ID']))
    intervals = natural_sort(set(df['Range']))
    return (intervals,)


@app.cell
def _(df, file_import, mo):
    mo.ui.table(
        df,
        label = f"##Smear Analysis\nFile: **{file_import.value[0].name}**",
        show_column_summaries=False,
        freeze_columns_left=["Sample ID"],
        show_data_types = False,
        selection = None,
        pagination = False,
        max_height=400
    )
    return


@app.cell
def _(intervals, mo):
    target_range = mo.ui.radio(intervals, value = intervals[-2], inline = True, label="Target range for the samples: ")
    target_pmol = mo.ui.slider(
        value = 15.0,
        start = 0.1,
        stop = 50.0,
        step = 0.1,
        include_input = True,
        full_width = True,
        label = "Target **picomoles** for final Libraries"
    )
    return target_pmol, target_range


@app.cell
def _(mo, target_pmol, target_range):
    mo.vstack([target_range, target_pmol])
    return


@app.cell
def _(df, mo, pd):
    input_header = mo.md("##Sample Concentrations\nUsing this interactive form, input the concentrations for each sample in **ng/µL** (nanograms per microliter).")

    def unique(sequence):
        seen = set()
        return [x for x in sequence if not (x in seen or seen.add(x))]
    samples_id = unique(df['Sample ID'])

    quants_df = mo.ui.data_editor(
        pd.DataFrame({
            'Sample ID' : list(samples_id),
            'concentration (ng/µL)': [0.0 for i in range(len(samples_id))]
        }),
        editable_columns= ['concentration (ng/µL)'],
        label = "Add sample concentrations here",
        pagination= False
    )

    mo.vstack([input_header,quants_df], gap = 0)
    return (quants_df,)


@app.cell
def _(
    df,
    intervals,
    mo,
    pd,
    process_sample,
    quants_df,
    target_pmol,
    target_range,
):
    mo.stop(any([i == 0 for i in quants_df.value['concentration (ng/µL)']]))

    n_rows = pd.DataFrame(df.groupby('Sample ID').size(),columns=['sample'])

    df_with_conc = df.merge(
        quants_df.value,
        left_on='Sample ID',
        right_on='Sample ID',
        how="left"
    )

    calc_table = df_with_conc.groupby('Well', sort = False).apply(lambda group: process_sample(group, target_range.value, intervals, target_pmol.value), include_groups=False)
    mo.ui.table(
        calc_table,
        pagination = False,
        selection = None,
        show_column_summaries = False,
        show_data_types = False,
        freeze_columns_left = ["Sample ID"],
        label = f"## Scaled Concentrations\n\nThis table scales the concentrations you input above with the proportion of the sample with the target `Range` **{target_range.value}**"
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
        full_width = True,
        label = "Volume to elute in (µL)"
    )
    elution_vol
    return (elution_vol,)


@app.cell
def _(calc_table, elution_vol, mo, pd, target_pmol):
    total_vol = calc_table['Volume to Pool'].sum()
    total_ng = calc_table['ng Primary Library'].sum()
    total_pM = target_pmol.value * len(calc_table.index)

    pool_ngul = 0 if total_vol == 0 else round(total_ng / total_vol,1)
    pool_uM = 0 if total_vol == 0 else total_pM / total_vol
    recovery = round(pool_uM * total_vol / elution_vol.value, 2)

    def style_cell(_rowId, _columnName, value):
        if _columnName == "Value":
            return {"fontWeight": "bold"}
        return {}

    mo.hstack(
        [mo.ui.table(
            pd.DataFrame({
                "Metric" : ["Total Volume of Pool (µL)", "Total ng Across Pools", "Pool ng/µL", "Total pM across intervals", "µM Per Pool across intervals", "µM Per Pool assuming 100% recovery"],
                "Value" : [round(total_vol,2), total_ng, round(pool_ngul,2) , total_pM, round(pool_uM,1), recovery]
            }),
            label = "## Final Pooling Metrics",
            show_data_types = False,
            style_cell = style_cell,
            pagination = False,
            selection = None
        ),
        ""],
        widths= [0, 1]
    )
    return


@app.cell
def _(file_import, mo, quants_df):
    mo.stop(not file_import.value or mo.stop(any([i == 0 for i in quants_df.value['concentration (ng/µL)']])))

    model_text = mo.md("----\n## Model and Scale\nUsing the data above, we can create a fitted model (non-linear least squares via log(x)) that will help scale the concentrations of the other libraries made alongside these that were not submitted for fragment analysis.")

    samples_import = mo.ui.file(
        kind="area",
        filetypes = [".csv", ".CSV", ".Csv"],
        label = "Import a CSV file of samples and their concentrations here."
    )

    headers = mo.ui.switch(value= True, label = "file has column headers")
    mo.vstack([model_text, samples_import, headers])
    return headers, samples_import


@app.cell
def _(file_import, headers, io, mo, pd, samples_import):
    mo.stop(not samples_import.value)

    try:
        libdf = pd.read_table(io.BytesIO(samples_import.value[0].contents), header= 0 if headers.value else None, engine='python', sep=None)
        libdf.columns = [i.lower() for i in libdf.columns]
    except Exception:
        is_err = True
        err_md = mo.md(f"""
        /// error| Invalid input file

        The `pandas` parser failed to read the input file. Please check that it conforms to a tabular comma delimited file. The first 200 bytes of the uploaded file:

        ///
        """)
        err = mo.vstack([err_md, mo.md(f"```\n{file_import.value[0].contents[:200]}\n```")])
        mo.stop(is_err, err)

    colnames = mo.ui.array([
        mo.ui.dropdown(options=list(libdf.columns), label=f"Sample Names"),
        mo.ui.dropdown(options=list(libdf.columns), label=f"Concentration"),
    ],
        label = "Column Names"
    )
    mo.vstack([
        mo.md("Please select the columns in your input file with the sample names and library concentrations."),
        mo.ui.table(libdf.head(3), selection = None, show_column_summaries=False, show_data_types=False, show_download = False),
        mo.hstack(colnames, justify="start")
    ])
    return colnames, libdf


@app.cell
def _(calc_table, colnames, libdf, mo, np, samples_import):
    mo.stop(
        not samples_import.value,
        mo.md(
            f"""
    /// admonition| Library concentration file required

    ///
            """
        )
    )
    #----------
    from sklearn.linear_model import LinearRegression

    x = np.log(calc_table["Sample ng/µL"].values).reshape(-1, 1)
    y = calc_table["Corrected ng/µL"].values

    fit = LinearRegression().fit(x, y)

    # Predict
    x_seq = libdf[colnames[-1].value].values
    y_pred = fit.predict(np.log(x_seq).reshape(-1, 1))
    return fit, x, x_seq, y, y_pred


@app.cell
def _(calc_table, fit, plt, x, x_seq, y, y_pred):
    rsq = f"R² = {fit.score(x, y):.3f}"
    plt.pyplot.scatter(calc_table["Sample ng/µL"], calc_table["Corrected ng/µL"], color="black", label="Frag Data")
    plt.pyplot.scatter(x_seq, y_pred, facecolors='none', edgecolors='dodgerblue', label="Corrected Libraries")
    plt.pyplot.xlabel("Sample ng/µL")
    plt.pyplot.ylabel("Corrected ng/µL")
    plt.pyplot.ylim(0, None)
    plt.pyplot.title(f"Model Fit ({rsq})")
    plt.pyplot.legend()
    plt.pyplot.gcf()
    return


@app.cell
def _(colnames, libdf, mo, pd, x_seq, y_pred):
    mo.ui.table(
        pd.DataFrame(
            {
                "Sample" : libdf[colnames[0].value].values,
                "Sample ng/µL": x_seq,
                "Frag-Corrected ng/µL" : y_pred.round(2)}
        ),
        selection=None,
        show_data_types=False,
        show_column_summaries=False
    )
    return


if __name__ == "__main__":
    app.run()
