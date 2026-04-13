import marimo

__generated_with = "0.21.1"
app = marimo.App(width="medium", app_title="Fragment Analysis Library Pooling")


@app.cell
def _():
    import io
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import re
    from sklearn.linear_model import LinearRegression
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

    sample_table = mo.md("""
    | Sample | Conc |
    |:-------|:-----|
    | sample_1 | 1.2 |
    | sample_2 | 0.74 |
    """)
    return (
        LinearRegression,
        example_table,
        io,
        mo,
        np,
        pd,
        plt,
        re,
        sample_table,
    )


@app.cell
def _(pd, re):
    def process_sample(group, target_identifier, all_intervals):
        target_row = group[group['Range'] == target_identifier].iloc[0]
        corrected_sum = group[group['Range'] != all_intervals[0]]['ng/µL'].sum()
        target_conc = target_row['ng/µL']
        target_size = target_row['Avg. Size']
        sample_id = target_row['Sample ID']
        corrected_smear = target_conc/corrected_sum
        quant = target_row['concentration (ng/µL)']
        nM = (quant * corrected_smear) / (target_size * 660) * 1000000
        return pd.Series(
            {
            'Sample ID': sample_id,
            'Avg.Size': target_size,
            '% of Total Conc.': round(corrected_smear * 100,1),
            'Window ng/µL': round(target_conc,3),
            'Sample ng/µL': quant,
            'Est. nM': round(nM, 2),
            'Corrected ng/µL': round(quant * corrected_smear,3)
        })

    def natural_sort(l): 
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(l, key=alphanum_key)

    return natural_sort, process_sample


@app.cell
def _(mo):
    file_import = mo.ui.file(
        kind="area",
        filetypes = [".csv", ".CSV", ".Csv"],
        label = "Import fragment analysis CSV file here"
    )

    samples_import = mo.ui.file(
        kind="button",
        filetypes = [".csv", ".CSV", ".Csv"],
        label = "Import"
    )

    sampleheaders = mo.ui.switch(value= True, label = "Sample file has column headers")
    return file_import, sampleheaders, samples_import


@app.cell
def _(file_import, mo, sampleheaders, samples_import):
    mo.sidebar(
        [
            mo.md('# Frag-Scaled Molarity\nThis worksheet scales the concentration/molarity of samples with the proportional to your target fragment interval as determined by fragment analysis.'),
            mo.vstack([mo.md("### Sample concentrations"),samples_import, sampleheaders]),
            file_import
        ],
        footer = mo.md('<img src="public/gih_logo.png" width="200" />\n\nMade with ❤️ for 🧬')
    )
    return


@app.cell
def _(mo, sample_table, samples_import):
    mo.md(f"""
    ## {'▷' if not samples_import.value else '✅'} Import Sample Concentrations
    Using the left sidebar, import a CSV file of samples and their concentrations. This file should  at minimum feature samples used in the fragment analysis.

    {mo.accordion({"🔍︎ Example Sample Concentration File" : sample_table})}
    """)
    return


@app.cell
def _(file_import, io, mo, pd, sampleheaders, samples_import):
    mo.stop(
        not samples_import.value,
        mo.md("/// admonition| Sample file required\n\nCannot proceed until the file is uploaded\n///")
    )

    try:
        sampdf = pd.read_csv(io.BytesIO(samples_import.value[0].contents), header= 0 if sampleheaders.value else None)
        sampdf.columns = [i.lower() for i in sampdf.columns]
    except Exception:
        is_err = True
        err_md = mo.md(f"""
        /// error| Invalid input file

        The `pandas` parser failed to read the input file. Please check that it conforms to a tabular comma delimited file. The first 200 bytes of the uploaded file:

        ///
        """)
        mo.stop(
            is_err,
            mo.vstack([err_md, mo.md(f"```\n{file_import.value[0].contents[:200]}\n```")])
        )
    return (sampdf,)


@app.cell
def _(example_table, file_import, mo):
    mo.vstack([
         mo.md(f"""## {'▷' if not file_import.value else '✅'} Import Smear Analysis\nUsing the left sidebar, import the CSV file that was provided to you by the BRC from the fragment analyzer. Then, choose which fragment size interval you are interested in performing a scaled concentration correction on."""),
        mo.accordion({
        "🔍︎ Example Fragment Analyzer CSV File": example_table,
        "⚠️ Notes and Considerations" : mo.md("**First row skipped**: Be aware that the first interval (sorted alphanumerically) of each sample, usually 10bp-100bp, is skipped in the calculations below.\n\n**Singletons skipped**: Rows with a `Sample ID` that only appears once are removed (e.g. a ladder, samples from a different run).\n\n**Consistency**: The target range is expected to be consistent across all samples.")})
    ])
    return


@app.cell
def _(file_import, io, mo, pd):
    mo.stop(
        not file_import.value,
        mo.md("/// admonition| Smear analysis file required\n\nCannot proceed until the file is uploaded\n///")
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
def _(df, mo, natural_sort, sampdf):
    intervals = natural_sort(set(df['Range']))
    target_range = mo.ui.radio(intervals, value = intervals[-2], inline = True)

    colnames = mo.ui.array([
        mo.ui.dropdown(options=list(sampdf.columns), label=f"Sample Names"),
        mo.ui.dropdown(options=list(sampdf.columns), label=f"Concentration"),
    ],
        label = "Column Names"
    )
    return colnames, intervals, target_range


@app.cell
def _(colnames, file_import, mo, samples_import, target_range):
    mo.stop(not file_import.value or not samples_import.value)

    mo.vstack([
        mo.md(f"### {'▷' if any([not i.value for i in colnames]) else '✅'} Configure Imports\nColumns in sample-concentration file:"),
        mo.hstack(colnames, justify="start"),
        mo.md("Target genomic size range for the samples:"),
        target_range
    ])
    return


@app.cell
def _(colnames, file_import, samples_import):
    imports_finished = not samples_import.value or any([not i.value for i in colnames]) or not file_import.value
    return (imports_finished,)


@app.cell
def _(colnames, df, imports_finished, mo, sampdf):
    mo.stop(imports_finished)
    quants_df = sampdf[[colnames[0].value, colnames[1].value]].rename(columns={colnames[0].value: 'Sample ID', colnames[1].value: 'concentration (ng/µL)'})

    lib_id = list(set(quants_df['Sample ID']))
    id_err = False

    try:
        sample_id = list(set(df['Sample ID']))
        if set(sample_id).intersection(set(lib_id)) != set(sample_id):
            raise ValueError
    except ValueError:
        id_err = True
        mo.stop(
            True,
            mo.md(f"""/// error| Sample mismatch\n
    The sample concentration file does not contain the same samples as the smear analysis file. Samples in smear analysis missing in concentration file:\n\n
    {'\n\n'.join(set(sample_id).difference(set(sample_id).intersection(set(lib_id))))}
    ///""")
        )
    return id_err, quants_df


@app.cell
def _(df, file_import, id_err, imports_finished, mo):
    mo.stop(imports_finished or id_err)
    mo.accordion({
        "View Smear Data": mo.ui.table(
        df,
        label = f"## Smear Analysis\nFile: **{file_import.value[0].name}**",
        show_column_summaries=False,
        freeze_columns_left=["Sample ID"],
        show_data_types = False,
        selection = None,
        pagination = False,
        max_height=200
        )
    }
    )
    return


@app.cell
def _(
    df,
    id_err,
    imports_finished,
    intervals,
    mo,
    pd,
    process_sample,
    quants_df,
    target_range,
):
    mo.stop(imports_finished or id_err)

    n_rows = pd.DataFrame(df.groupby('Sample ID').size(),columns=['sample'])

    df_with_conc = df.merge(
        quants_df,
        left_on='Sample ID',
        right_on='Sample ID',
        how="left"
    )

    calc_table = df_with_conc.groupby('Well', sort = False).apply(lambda group: process_sample(group, target_range.value, intervals), include_groups=False)
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
def _(id_err, imports_finished, mo):
    mo.stop(imports_finished or id_err)

    mo.md("----\n## Model and Scale\nUsing the data above, we can create a fitted model (non-linear least squares via log(x)) that will help scale the concentrations of the other libraries made alongside these that were not submitted for fragment analysis.")
    return


@app.cell
def _(x_seq):
    min(x_seq)
    return


@app.cell
def _(
    LinearRegression,
    calc_table,
    id_err,
    imports_finished,
    mo,
    np,
    quants_df,
):
    mo.stop(imports_finished or id_err)

    predict_table = calc_table[calc_table["Sample ng/µL"] >= 0.2]

    x = np.log(predict_table["Sample ng/µL"].values).reshape(-1, 1)
    linear_x = predict_table["Sample ng/µL"].values.reshape(-1, 1)
    y = predict_table["Corrected ng/µL"].values
    nMy = predict_table["Est. nM"].values
    Sizey = predict_table["Avg.Size"].values

    fit = LinearRegression().fit(x, y)
    Molarityfit = LinearRegression().fit(x, nMy)
    Sizefit = LinearRegression().fit(x, Sizey)

    # Predict
    x_seq = quants_df['concentration (ng/µL)'].values

    concidentity_pred = x.flatten()
    concss_res = np.sum((y - concidentity_pred) ** 2)
    concss_tot = np.sum((y - np.mean(y)) ** 2)
    concidentity_r2 = 1 - concss_res / concss_tot


    y_pred = fit.predict(np.log(x_seq).reshape(-1, 1))
    rsq = f"Log R² = {fit.score(x, y):.3f} | 1:1 Linear R² = {concidentity_r2:.3f}"

    identity_pred = linear_x.flatten()
    ss_res = np.sum((nMy - identity_pred) ** 2)
    ss_tot = np.sum((nMy - np.mean(nMy)) ** 2)
    identity_r2 = 1 - ss_res / ss_tot

    nMy_pred = Molarityfit.predict(np.log(x_seq).reshape(-1, 1))
    nM_rsq = f"Log R² = {Molarityfit.score(x, nMy):.3f} | 1:1 Linear R² = {identity_r2:.3f}"

    size_pred = Sizefit.predict(np.log(x_seq).reshape(-1, 1))
    return nM_rsq, nMy_pred, rsq, size_pred, x_seq, y_pred


@app.cell
def _(
    calc_table,
    id_err,
    imports_finished,
    mo,
    nM_rsq,
    nMy_pred,
    plt,
    rsq,
    x_seq,
    y_pred,
):
    mo.stop(imports_finished or id_err)

    fig, axes = plt.subplots(2,1,sharex=True, figsize=(10, 7))

    # Plot on the first axes
    axes[0].set_title(f"Concentration ({rsq})")
    axes[0].scatter(calc_table["Sample ng/µL"], calc_table["Corrected ng/µL"], color="black", label="Frag Data")
    axes[0].scatter(x_seq, y_pred, facecolors='none', edgecolors='dodgerblue', label="Modeled Libraries")
    axes[0].axline((0,0), slope=1, color = 'indianred', dashes = (1, 2, 1, 2))
    axes[0].set_xlabel("Sample ng/µL")
    axes[0].set_ylabel("Smear-Corrected ng/µL")
    axes[0].set_ylim(0, None)
    axes[0].legend()

    # Plot on the second axes
    axes[1].set_title(f"Molarity ({nM_rsq})")
    axes[1].scatter(calc_table["Sample ng/µL"], calc_table["Est. nM"], color="black", label="Frag Data")
    axes[1].scatter(x_seq, nMy_pred, facecolors='none', edgecolors='darkseagreen', label="Modeled Libraries")
    axes[1].axline((0,0), slope=1, color = 'indianred', dashes = (1, 2, 1, 2))
    axes[1].set_xlabel("Sample ng/µL")
    axes[1].set_ylabel("Estimated nM")
    axes[1].set_ylim(0, None)
    axes[1].legend()

    # 4. Use tight_layout to prevent label overlap
    plt.tight_layout()
    plt.gcf()
    return


@app.cell
def _(nMy_pred, pd, quants_df, size_pred, x_seq, y_pred):
    def warn_cell(row_id, column_name, value):
        # row_id is a string index; look up the target column value for this row
        try:
            row_val = frag_scaled_concs.iloc[int(row_id)]["Estimated nM"]
            if isinstance(row_val, (int, float)):
                if row_val <= 0:
                    return {"backgroundColor": "lightcoral", "color": "white"}
                elif row_val < 0.3:
                    return {"backgroundColor": "#f3a155", "color": "white"}
        except (ValueError, KeyError, IndexError):
            pass
        return {}

    frag_scaled_concs = pd.DataFrame(
            {
                "Sample" : quants_df['Sample ID'].values,
                "Quant ng/µL": x_seq,
                "Frag-Corrected ng/µL" : y_pred.round(2),
                "Average Fragment Size": size_pred.round(),
                "Estimated nM" : nMy_pred.round(2)
            }
        )
    return frag_scaled_concs, warn_cell


@app.cell
def _(frag_scaled_concs, mo, warn_cell):
    dropouts = sum([i <= 0 for i in frag_scaled_concs["Estimated nM"].values])
    likelydropouts = sum([(i > 0 and i < 0.3) for i in frag_scaled_concs["Estimated nM"].values])

    mo.vstack([
        mo.md(f"Obvious dropouts: **{dropouts}** (red) | Likely dropouts: **{likelydropouts}** (orange)"),
        mo.ui.table(
            frag_scaled_concs,
            page_size = 24,
            show_data_types=False,
            show_column_summaries=False,
            style_cell=warn_cell
        )
    ])
    return


if __name__ == "__main__":
    app.run()
