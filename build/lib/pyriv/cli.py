# -*- coding: utf-8 -*-

"""Console script for pyriv."""

import click
from graph_prep import GraphBuilder
from river_graph import RiverGraph


@click.command()
def main(args=None):
    """Console script for pyriv."""
    click.echo("this will be the command line interface script for pyriv. this will be the entry point to the program")
    click.echo("See click documentation at http://click.pocoo.org/")


if __name__ == "__main__":
    main()
