import pandas as pd
import re

df = pd.read_parquet("wikitext-2/train.parquet")

words = []
for line in df["text"]:
    for token in line.split():
        if re.fullmatch(r"[A-Za-z]+", token):  # only pure alphabetic words
            words.append(token.lower())        # lowercase

# Save to one-word-per-line text file
with open("wikitext-2-raw.txt", "w") as f:
    f.write("\n".join(words))
