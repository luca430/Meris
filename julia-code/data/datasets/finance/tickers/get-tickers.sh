#!/bin/bash
# Get lists of tickers for relevant markets
set -e

NYSEURL=https://datahub.io/core/nyse-other-listings/_r/-/data/nyse-listed.csv
curl -L $NYSEURL -o "nyse-tickers.csv"

NASDAQURL=https://datahub.io/core/nasdaq-listings/_r/-/data/nasdaq-listed.csv
curl -L $NASDAQURL -o "nasdaq-tickers.csv"

# note: EURONEXT tickers may be obtained at
#       https://live.euronext.com/en/products/equities/list
