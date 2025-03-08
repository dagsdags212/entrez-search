import marimo

__generated_with = "0.11.17"
app = marimo.App(width="medium")


@app.cell
def _():
    from pathlib import Path
    import asyncio

    import marimo as mo
    import polars as pl
    import altair as alt
    from great_tables import GT
    from Bio import Entrez
    from Bio import SeqIO

    Entrez.email = "jegsamson.dev@gmail.com"
    Entrez.tool = "biopython"
    return Entrez, GT, Path, SeqIO, alt, asyncio, mo, pl


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

    retmax = mo.ui.number(
        start=5, stop=10000, step=10,
        value=50,
        debounce=True,
        label="Limit records to:",
        full_width=True,
    )
    return db_dropdown, db_list, query, retmax


@app.cell
def _(db_dropdown, mo, query, retmax):
    mo.vstack([
        mo.md("# Entrez Query Builder").center(),
        mo.hstack([
            db_dropdown, 
            query,
        ], align="stretch", widths=[1, 10]),
        retmax,
        mo.md(f"From the **{db_dropdown.value}** database, I want to search for the term").center(),
        mo.md(f"## {mo.icon('svg-spinners:bouncing-ball')} {query.value} {mo.icon('svg-spinners:bouncing-ball')}").center(),
        mo.md(f"while limiting the results to **{retmax.value}** records.").center()
    ])
    return


@app.cell
def _(Entrez, Path, SeqIO, mo, pl):
    def fetch_records(query: str, db: str, retmax: int) -> pl.DataFrame:
        esearch_handle = Entrez.esearch(db=db, term=query, retmax=retmax)
        id_list = Entrez.read(esearch_handle).get("IdList")
        results = Entrez.read(Entrez.esummary(db=db, id=id_list))
        return pl.DataFrame(results).select(
            pl.col("Id").cast(pl.String),
            pl.col("TaxId").cast(pl.String),
            pl.col(["AccessionVersion", "Caption", "Title", "Length"]),
            pl.col("CreateDate").str.to_date("%Y/%m/%d")
        )

    def download_genbank_record(id: int | str, progress_bar: mo.status.progress_bar) -> None:
        """Downloads a single GenBank record from a given identifier."""
        handle = Entrez.efetch(db="nuccore", id=id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "gb")
    
        outdir = Path("data")
        if not outdir.exists():
            outdir.mkdir(exist_ok=True)

        outfile = outdir / f"{record.id}.gbk"
        if outfile.exists():
            progress_bar.update(title="Skipping file...", subtitle=f"FASTA file for {id} already exists")
            return

        SeqIO.write(record, outfile, "gb")

    def download_fasta_file(id: int | str, progress_bar: mo.status.progress_bar) -> None:
        """Downloads a single FASTA file from a given identifier."""
        handle = Entrez.efetch(db="nuccore", id=id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")

        outdir = Path("data")
        if not outdir.exists():
            outdir.mkdir(exist_ok=True)

        outfile = outdir / f"{record.id}.fasta"
        if outfile.exists():
            progress_bar.update(title="Skipping file...", subtitle=f"FASTA file for {id} already exists")
            return

        SeqIO.write(record, outfile, "fasta")
    

    def batch_download_gb(id_list: list) -> None:
        """Retrieve a set of GenBank records using Entrez esearch."""
        with mo.status.progress_bar(
            total=len(id_list),
            title="Downloading GenBank files",
            subtitle="Please wait",
            show_rate=True, show_eta=True,
            completion_title=f"Done!"
        ) as bar: 

            for id in id_list:
                bar.update(title="Downloading files", subtitle=f"Fetching GBK file for {id}")
                download_genbank_record(id, bar)
            
    def batch_download_fasta(id_list: list) -> None:
        """Retrieve a set of GenBank records using Entrez esearch."""
        with mo.status.progress_bar(
            total=len(id_list),
            title="Downloading FASTA files",
            subtitle="Please wait",
            show_rate=True, show_eta=True,
            completion_title=f"Done!"
        ) as bar: 

            for id in id_list:
                bar.update(title="Downloading files", subtitle=f"Fetching FASTA file for {id}")
                download_fasta_file(id, bar)
    return (
        batch_download_fasta,
        batch_download_gb,
        download_fasta_file,
        download_genbank_record,
        fetch_records,
    )


@app.cell
def _(db_dropdown, fetch_records, query, retmax):
    df = fetch_records(query.value, db_dropdown.value, retmax.value)
    return (df,)


@app.cell
def _(df, mo):
    # Filter out outliers.
    lengths = df.select("Length")
    dates = df.select("CreateDate").to_series()

    length_filter = mo.ui.range_slider(
        label="By length",
        start=lengths.quantile(0.25).item(),
        stop=lengths.quantile(0.75).item(),
        full_width=True,
        show_value=True,
    )

    date_filter = mo.ui.date_range.from_series(
        dates,
        label="By date",
        full_width=True,
    )
    return date_filter, dates, length_filter, lengths


@app.cell
def _(alt, date_filter, df, length_filter, pl):
    sql_query = f"""
        SELECT *
        FROM self
        WHERE 
            Length BETWEEN {length_filter.value[0]} AND {length_filter.value[1]}
                OR 
            CreateDate BETWEEN {date_filter.value[0]} AND {date_filter.value[1]};
    """

    filtered_df = df.sql(sql_query)

    min_length = filtered_df.select("Length").min().item()
    max_length = filtered_df.select("Length").max().item()
    n_unqiue_taxid = len(filtered_df.select("TaxId").unique())

    taxid_counts = filtered_df.select("TaxId").group_by("TaxId").agg(pl.len().alias("Count"))

    taxid_counts_chart = alt.Chart(taxid_counts).mark_arc().encode(
        theta="Count", color="TaxId"
    )
    return (
        filtered_df,
        max_length,
        min_length,
        n_unqiue_taxid,
        sql_query,
        taxid_counts,
        taxid_counts_chart,
    )


@app.cell
def _(
    date_filter,
    filtered_df,
    length_filter,
    max_length,
    min_length,
    mo,
    taxid_counts_chart,
):
    filters = mo.vstack([
        mo.md("## Filters"),
        mo.hstack([date_filter, length_filter], widths=[1, 4], gap=3,align="center"),
    ])


    mo.vstack([
        filters,
        mo.md(f"## {mo.icon('fluent-mdl2:poll-results')}   Results"),
        mo.hstack([
            mo.vstack([    
                mo.stat(value=len(filtered_df.select("TaxId").unique()), label="Unique Taxa", bordered=True),
                mo.stat(value=min_length, label="Minimum Length", bordered=True),
                mo.stat(value=max_length, label="Maximum Length", bordered=True),
            ]),
            mo.ui.altair_chart(taxid_counts_chart)
        ], widths=[1, 4], align="center"), 
    ])
    return (filters,)


@app.cell
def _(filtered_df, mo):
    table = mo.ui.table(data=filtered_df, pagination=True)
    return (table,)


@app.cell
def _(fa_download_btn, gb_download_btn, mo, table):
    mo.vstack([
        mo.md(f"## {mo.icon('fluent-color:table-24')} Records Table"),
        table,
        mo.hstack([
            mo.icon('line-md:downloading-loop'),
            fa_download_btn, 
            gb_download_btn
        ], justify="start", align="center") 
    ])
    return


@app.cell
def _(batch_download_fasta, batch_download_gb, mo, table):
    gb_download_btn = mo.ui.button(
        value=table.value["Id"].to_list(), 
        label="Download GBK",
        kind="success",
        tooltip="Download GenBank Records",
        on_click=batch_download_gb,
    )

    fa_download_btn = mo.ui.button(
        value=table.value["Id"].to_list(),
        label="Download FASTA",
        kind="info",
        tooltip="Download FASTA files",
        on_click=batch_download_fasta,
    )
    return fa_download_btn, gb_download_btn


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
