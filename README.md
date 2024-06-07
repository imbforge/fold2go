# fold2go

## Description

This is a Nextflow pipeline to run AlphaFold on IMB infrastructure. Multiple sequence alignments have been factored out to avoid blocking GPU resources with CPU workloads.
It implements a Shiny (Python) application that allows to track pipeline progress and results interactively.

## Overview

```mermaid
flowchart LR
    subgraph fold2go
        subgraph "[GPU]"
        v17([ALPHAFOLD])
        end
        subgraph "[CPU]"
        v11([MSA])
        v19([PYMOL])
        v7([SHINY])
        end
    v11 --> v17
    v17 --> v19
    v19 --> v7
    v17 --> v7
    end
```

## Graphical User Interface

This pipeline can be used via an in-development Jupyterhub-based graphical frontend. If you want to try this out, head over to [imb-alphafold](https://gitlab.rlp.net/imbforge/imb-alphafold)
