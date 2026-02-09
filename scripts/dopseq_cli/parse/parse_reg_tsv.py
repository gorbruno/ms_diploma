#!/usr/bin/env python3
from pathlib import Path
import click
import pandas as pd
import logging
import sys

LOG_FILE = ".log"

def setup_logging(verbose: bool):
    handlers = [logging.FileHandler(LOG_FILE, mode="a")]
    if verbose:
        handlers.append(logging.StreamHandler())

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=handlers
    )
    return logging.getLogger(__name__)


@click.command()
@click.option(
    "--input", "-i",
    "input_dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Directory with *.reg.tsv files."
)
@click.option(
    "--outdir", "-o",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True),
    help="Directory to store all output files."
)
@click.option(
    "--pattern", "-p",
    default="*.reg.tsv",
    help="File search pattern (default: *.reg.tsv)."
)
@click.option(
    "--reads-min", "-r",
    type=int,
    default=None,
    help="Minimum reg_reads value. If omitted — filter is not applied."
)
@click.option(
    "--prefix-contigs",
    default="C",
    help="Prefix to classify chromosome as contigs (default: C)."
)
@click.option(
    "--prefix-chroms",
    default="O",
    help="Prefix to classify chromosome as chroms (default: O)."
)
@click.option(
    "--excel",
    is_flag=True,
    help="Also save processed tables to Excel format."
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    help="Enable verbose logging (logs also printed to console)."
)
async def main(input_dir, outdir, pattern, reads_min, prefix_contigs, prefix_chroms, excel, verbose):
    logger = setup_logging(verbose)

    # Log the command itself
    command_str = " ".join(sys.argv)
    logger.info(f"Executed command: {command_str}")

    input_dir = Path(input_dir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Input directory: {input_dir}, Output directory: {outdir}")

    excel_dir = outdir / "excel"
    if excel:
        excel_dir.mkdir(exist_ok=True)
        logger.info(f"Excel output enabled, directory: {excel_dir}")

    files = sorted(list(input_dir.glob(pattern)))
    if not files:
        msg = f"No files found matching pattern '{pattern}' in {input_dir}"
        click.echo(msg)
        logger.warning(msg)
        return

    click.echo(f"Found files: {len(files)}")
    logger.info(f"Found {len(files)} files matching pattern '{pattern}'")

    # Storage for summary Excel
    summary_contigs = {}
    summary_chroms = {}

    for file in files:
        file_label = file.name.replace(".reg.tsv", "")
        click.echo(f"\nProcessing: {file.name}")
        logger.info(f"Processing file: {file}")

        df = pd.read_csv(file, sep="\t")
        df = df.sort_values("p_value", ascending=True)

        if reads_min is not None:
            df = df[df["reg_reads"] > reads_min]
            logger.info(f"Filtered by reg_reads > {reads_min}, remaining rows: {len(df)}")

        contigs_df = df[df["chrom"].astype(str).str.startswith(prefix_contigs)]
        chroms_df = df[df["chrom"].astype(str).str.startswith(prefix_chroms)]

        out_contigs = outdir / f"{file_label}.contigs.tsv"
        out_chroms = outdir / f"{file_label}.chroms.tsv"

        contigs_df.to_csv(out_contigs, sep="\t", index=False)
        chroms_df.to_csv(out_chroms, sep="\t", index=False)
        click.echo(f"→ Saved: {out_contigs.name}, {out_chroms.name}")
        logger.info(f"Saved TSVs: {out_contigs}, {out_chroms}")

        summary_contigs[file_label] = contigs_df
        summary_chroms[file_label] = chroms_df

        if excel:
            excel_contigs = excel_dir / f"{file_label}.contigs.xlsx"
            excel_chroms = excel_dir / f"{file_label}.chroms.xlsx"
            contigs_df.to_excel(excel_contigs, index=False)
            chroms_df.to_excel(excel_chroms, index=False)
            click.echo(f"→ Excel: {excel_contigs.name}, {excel_chroms.name}")
            logger.info(f"Saved Excel files: {excel_contigs}, {excel_chroms}")

    if excel:
        summary_contigs_path = excel_dir / "summary_contigs.xlsx"
        summary_chroms_path = excel_dir / "summary_chroms.xlsx"

        with pd.ExcelWriter(summary_contigs_path) as writer:
            for sheet_name, sdf in summary_contigs.items():
                sdf.to_excel(writer, sheet_name=sheet_name[:31], index=False)

        with pd.ExcelWriter(summary_chroms_path) as writer:
            for sheet_name, sdf in summary_chroms.items():
                sdf.to_excel(writer, sheet_name=sheet_name[:31], index=False)

        click.echo("\n→ Summary Excel files created:")
        click.echo(f"   - {summary_contigs_path.name}")
        click.echo(f"   - {summary_chroms_path.name}")
        logger.info(f"Created summary Excel files: {summary_contigs_path}, {summary_chroms_path}")

    click.echo("\nDone!")
    logger.info("Processing completed successfully.")


if __name__ == "__main__":
    main()