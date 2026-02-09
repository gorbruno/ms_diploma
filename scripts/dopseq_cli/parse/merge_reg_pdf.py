#!/usr/bin/env python3
from pathlib import Path
import click
from PyPDF2 import PdfMerger
import logging
import sys

LOG_FILE = ".log"

def setup_logging(verbose: bool):
    handlers = [
        logging.FileHandler(LOG_FILE, mode="a")
    ]

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
    help="Directory containing *.reg.pdf files."
)
@click.option(
    "--out", "-o",
    type=click.Path(dir_okay=False, file_okay=True),
    default="merged.pdf",
    help="Output merged PDF filename (default: merged.pdf)."
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    help="Enable verbose logging (logs also printed to console)."
)
async def main(input_dir, out, verbose):
    logger = setup_logging(verbose)

    # Log the command itself
    command_str = " ".join(sys.argv)
    logger.info(f"Executed command: {command_str}")

    input_dir = Path(input_dir)
    out = Path(out)

    logger.info(f"Starting merge. Input directory: {input_dir}, Output: {out}")

    pdfs = sorted(input_dir.glob("*.reg.pdf"))

    if not pdfs:
        msg = "No *.reg.pdf files found. Nothing to merge."
        click.echo(msg)
        logger.warning(msg)
        return

    click.echo(f"Found {len(pdfs)} PDF files.")
    logger.info(f"Found {len(pdfs)} PDF files to merge.")

    merger = PdfMerger()

    for pdf in pdfs:
        click.echo(f"Adding: {pdf.name}")
        logger.info(f"Adding file: {pdf}")
        with pdf.open("rb") as f:
            merger.append(f)

    with out.open("wb") as fout:
        merger.write(fout)

    msg = f"Saved merged PDF: {out}"
    click.echo(f"\n→ {msg}")
    click.echo("Done!")
    logger.info(msg)
    logger.info("Finished successfully.")


if __name__ == "__main__":
    main()