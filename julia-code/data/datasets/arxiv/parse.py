import os
from pathlib import Path
from nltk.stem import WordNetLemmatizer
wnl = WordNetLemmatizer()

ROOT = Path(os.path.join(os.path.dirname(__file__), 'raw-text/'))

#~ Loop through all available files and stem
for path in ROOT.rglob('*'):
    if path.is_file():
        print(f"Lemmatizing {path.name}...", end="\r", flush=True)
        NEWROOT = Path(str(path.parent).replace('raw-text', 'lemmatized-text'))
        NEWROOT.mkdir(parents=True, exist_ok=True)
        NEWPATH = NEWROOT / path.name

        # Parse the actual file
        #~Load
        words = Path(path).read_text().splitlines()
        #~Parse
        singles = [wnl.lemmatize(word) for word in words]
        #~Save
        Path(NEWPATH).write_text("\n".join(singles))
            
