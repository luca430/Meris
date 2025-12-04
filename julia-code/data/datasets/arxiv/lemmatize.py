import os
import numpy as np
import nltk
from nltk import pos_tag, word_tokenize
from nltk.corpus import wordnet
from nltk.stem import WordNetLemmatizer

nltk.download('punkt')
nltk.download('averaged_perceptron_tagger')

lemm = WordNetLemmatizer()

def wn_pos(treebank_tag):
    if treebank_tag.startswith('J'):
        return wordnet.ADJ
    elif treebank_tag.startswith('V'):
        return wordnet.VERB
    elif treebank_tag.startswith('N'):
        return wordnet.NOUN
    elif treebank_tag.startswith('R'):
        return wordnet.ADV
    else:
        return None

LEVEL_MAP = {
    "noun": wordnet.NOUN,
    "verb": wordnet.VERB,
    "adj":  wordnet.ADJ,
    "adv":  wordnet.ADV,
}

def lemmatize(sentence, level="noun"):
    requested = level.split("-")
    allowed_tags = {LEVEL_MAP[l] for l in requested if l in LEVEL_MAP}
    tokens = word_tokenize(sentence)
    tagged = pos_tag(tokens)
    out = []
    for w, pos in tagged:
        tag = wn_pos(pos)
        if tag in allowed_tags:
            out.append(lemm.lemmatize(w, tag))
        else:
            out.append(w)
    return out

# Paths
raw_folder = "clean_text"
clean_folder = "lemmatized_noun_adj_data"

for root, dirs, files in os.walk(raw_folder):
    for f in files:
        # Full path to input file
        infile = os.path.join(root, f)
        
        # Read file (assuming text files)
        try:
            txt = np.loadtxt(infile, dtype=str, encoding='utf-8')
        except:
            print(f"Skipping {infile} (cannot read as text)")
            continue
        
        # Join text and lemmatize
        try:
            iter(txt)
        except:
            print(f"Skipping {infile} (loaded file is not iterable)")
            continue
            
        text = " ".join(txt)
        lemmas = lemmatize(text, level="noun-adj")  # choose levels you want. OPTIONS: "noun-verb-adj-adv"
        
        # Compute output path
        rel_path = os.path.relpath(root, raw_folder)
        outdir = os.path.join(clean_folder, rel_path)
        os.makedirs(outdir, exist_ok=True)
        outfile = os.path.join(outdir, f)
        
        # Save lemmatized text, one lemma per line
        with open(outfile, 'w', encoding='utf-8') as out:
            if len(lemmas) != 0:
                for lemma in lemmas:
                    out.write(lemma + "\n")
        
        print(f"Processed {infile} → {outfile}")
