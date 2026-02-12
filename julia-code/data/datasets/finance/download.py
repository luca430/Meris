#!/usr/bin/env python3
"""
Download financial data from Yahoo finance using `yfinance`

Requires list of tickers at `tickers/tickers.csv`
Currently defaults to a list of NASDAQ tickers, from:
    https://datahub.io/core/nasdaq-listings#nasdaq-listed
"""
import os
import pandas as pd
import yfinance as yf

#~ Suppress some warnings related to dropping timezone information
#  When we assign "sample-id", we do not care about timezones, so we can safely drop it.
#  In addition, as we only use NASDAQ stocks [by default], timezones are the same for all.
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='pandas')

"Get NASDAQ tickers"
def get_nasdaq_tickers(filename: str = "tickers/nasdaq-tickers.csv"):
    tickers = pd.read_csv(filename)["ACT Symbol"].dropna().astype(str).tolist()
    return tickers

"Get NYSE tickers"
def get_nyse_tickers(filename: str = "tickers/nyse-tickers.csv"):
    tickers = pd.read_csv(filename)["ACT Symbol"].dropna().astype(str).tolist()
    return tickers

"""
Euronext tickers are a bit special as they need to be augmented with a suffix that specifies
the market. For example, Paris exchange gets suffix ".PA".
"""
def get_euronext_tickers(filename: str = "tickers/euronext-tickers.csv"):
    euronext_suffix = {"Paris": ".PA", "Milan": ".MI", "Amsterdam": ".AS"}
    
    with open(filename) as f:
        #~ Skip some lines that do not contain the separator
        lines = [line for line in f if ";" in line]
        #~ Load
        #  note: We have to provide `decimal=","` as European standards have decimal separators
        #        for their currencies and the ticker list contains the latest stock prices.
        df = pd.read_csv(pd.io.common.StringIO("".join(lines)), sep=";", decimal=",")

    #/ Append proper suffix for `yfinance`
    def detect_suffix(market):
        for key, suffix in euronext_suffix.items():
            if key in market:
                return suffix
        return None
    df["suffix"] = df["Market"].apply(detect_suffix)
    df = df[df["suffix"].notna()]
    tickers = (df["Symbol"] + df["suffix"]).tolist()
    
    #~ Return tickers and name
    return tickers

def get_volumes(tickers: str):
    #~ Allocate
    stockdataframes = {"h": [], "D": [], "W": [], "M": []}

    for ticker in tickers:
        print(f"Getting volumes for {ticker}...", end="\r")
        stockhistory = (
            yf.Ticker(ticker)
            .history(interval="1h", start="2025-01-01", end="2026-01-01", auto_adjust=False)
        )

        # Skip if empty
        if stockhistory.empty:
            continue
        
        # Resample; hourly, daily, weekly, montly
        aggregations = {
            "h": stockhistory["Volume"],
            "D": stockhistory["Volume"].resample("D").sum(),
            "W": stockhistory["Volume"].resample("W").sum(),
            "M": stockhistory["Volume"].resample("ME").sum(),
        }

        for period, volumes in aggregations.items():
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="Converting to PeriodArray/Index representation*.",
                    category=UserWarning,
                )
                df = pd.DataFrame({
                    "sample-id": volumes.index.tz_localize(None).to_period(period),
                    "ticker": ticker,
                    "total_volume": volumes.values,
                })
                stockdataframes[period].append(df)

    hdf = pd.concat(stockdataframes["h"], ignore_index=True)
    ddf = pd.concat(stockdataframes["D"], ignore_index=True)
    wdf = pd.concat(stockdataframes["W"], ignore_index=True)
    mdf = pd.concat(stockdataframes["M"], ignore_index=True)
    stockdataframes = {"hourly": hdf, "daily": ddf, "weekly": wdf, "monthly": mdf}
            
    return stockdataframes

def main(save: bool = True, out: str = "raw-data/") -> pd.DataFrame:
    #~ Note: to get tickers, see `tickers/get-tickers.sh`
    markets = ["euronext", "nasdaq", "nyse"]
    tickerfunctions = [get_euronext_tickers, get_nasdaq_tickers, get_nyse_tickers]
    for (i, market) in enumerate(markets):
        tickers = tickerfunctions[i](f"tickers/{market}-tickers.csv")
        #~ Extract volumes
        volumedf = get_volumes(tickers=tickers)
        #~ Save [when desired]
        if save:
            os.makedirs(os.path.dirname(out), exist_ok=True)
            for period in volumedf.keys():
                _df = volumedf[period]
                _df.to_csv(out + f"{market}-{period}-volumes.csv", sep=",", index=False)

if __name__ == "__main__":
    main()
