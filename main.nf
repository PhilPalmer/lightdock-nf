#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}