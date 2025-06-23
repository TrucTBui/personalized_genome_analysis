#!/bin/bash

# This script is used to lift over the coordinates of a pangenome analysis from one genome assembly to another using the UCSC liftOver tool.
CHAIN_FILE="~/tools/LiftOver/hg38ToHg19.over.chain.gz"
LIFTOVER="~/tools/LiftOver/liftOver"
INPUT_FILE=$1
LIFTOVER_OUTPUT=$2
UNMAPPED_OUTPUT=$3

# Check if the input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Input file not found!"
    exit 1
fi

# Check if the chain file exists
if [ ! -f "$CHAIN_FILE" ]; then
    echo "Chain file not found!"
    exit 1
fi

# Check if the liftOver tool exists
if [ ! -f "$LIFTOVER" ]; then
    echo "liftOver tool not found!"
    exit 1
fi

mkdir -p "$(dirname "$LIFTOVER_OUTPUT")"
mkdir -p "$(dirname "$UNMAPPED_OUTPUT")"

# Run liftOver
$LIFTOVER "$INPUT_FILE" "$CHAIN_FILE" "$LIFTOVER_OUTPUT" "$UNMAPPED_OUTPUT" -bedPlus=3

# Check if liftOver was successful
if [ $? -ne 0 ]; then
    echo "liftOver failed for '$INPUT_FILE'!"
    exit 1
fi

#Delete the unmapped output file if it is empty
if [ -s "$UNMAPPED_OUTPUT" ]; then
else
    echo "Unmapped output file "$UNMAPPED_OUTPUT" is empty. Deleting it."
    rm "$UNMAPPED_OUTPUT"
fi

