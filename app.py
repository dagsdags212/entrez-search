import marimo

__generated_with = "0.11.17"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    import altair as alt
    from Bio import Entrez

    Entrez.email = "jegsamson.dev@gmail.com"
    Entrez.tool = "biopython"
    return Entrez, alt, mo, pl


@app.cell
def _(mo):
    mo.md("""
    # Entrez Search
    """).center()
    return


@app.cell
def _(Entrez, mo):
    db_list = Entrez.read(Entrez.einfo())["DbList"]
    db_dropdown = mo.ui.dropdown(
       options=db_list,
        value="nuccore",
        label="",
        searchable=True,
        full_width=True,
    )

    query = mo.ui.text(
        value="Ebola virus",
        label="",
        placeholder="Query string",
        kind="text",
        full_width=True
    )

    retmax = mo.ui.slider(
        start=5, stop=10000, step=10,
        debounce=True,
        orientation="horizontal",
        show_value=True,
        label="Limit records to:",
        full_width=True,
    )
    return db_dropdown, db_list, query, retmax


@app.cell
def _(db_dropdown, mo, query, retmax):
    mo.vstack([
        mo.md("## Query Builder").center(),
        mo.hstack([
            db_dropdown, 
            query,
        ], align="stretch", widths=[1, 10]),
        retmax,
        mo.md(f"From the **{db_dropdown.value}** database, I want to search for the term ***{query.value}*** while limiting the results to **{retmax.value}** records.").center().callout(kind="info")
    ])
    return


@app.cell
def _(Entrez, pl):
    def fetch_records(query: str, db: str, retmax: int) -> pl.DataFrame:
        esearch_handle = Entrez.esearch(db=db, term=query, retmax=retmax)
        id_list = Entrez.read(esearch_handle).get("IdList")
        results = Entrez.read(Entrez.esummary(db=db, id=id_list))
        return pl.DataFrame(results).select(
            pl.col("Id").cast(pl.String),
            pl.col("TaxId").cast(pl.String),
            pl.col(["AccessionVersion", "Caption", "Title", 
             "CreateDate", "Length"])
        )

    return (fetch_records,)


@app.cell
def _(db_dropdown, fetch_records, query, retmax):
    df = fetch_records(query.value, db_dropdown.value, retmax.value)
    return (df,)


@app.cell
def _(df, pl):
    taxid_counts = df.select("TaxId").group_by("TaxId").agg(pl.len().alias("Count"))
    return (taxid_counts,)


@app.cell
def _(alt, df, mo, taxid_counts):
    _taxid_counts_chart = alt.Chart(taxid_counts).mark_arc().encode(
        theta="Count", color="TaxId"
    )

    mo.vstack([
        mo.md("## Records Table"),
        df, 
        mo.hstack([
            mo.ui.altair_chart(_taxid_counts_chart)
        ])
    ])

    return


@app.cell
def _(df):
    df.select()
    return


if __name__ == "__main__":
    app.run()
