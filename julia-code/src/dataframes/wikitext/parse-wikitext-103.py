import pandas as pd
import re

df1 = pd.read_parquet("wikitext-103/train-00000-of-00002.parquet")
df2 = pd.read_parquet("wikitext-103/train-00001-of-00002.parquet")

words = []
df = pd.concat([df1, df2], ignore_index=True)
for line in df["text"]:
    for token in line.split():
        if re.fullmatch(r"[A-Za-z]+", token):  # only pure alphabetic words
            words.append(token.lower())        # lowercase

# Save to one-word-per-line text file
with open("wikitext-103-raw.txt", "w") as f:
    f.write("\n".join(words))
