#~ Read a Chinese UTF8 encoded `.txt` and output it in the proper encoding
#  This will make parsing it down the line much easier
with open('story-of-the-stone.txt', encoding='utf-8', errors='ignore') as f, \
     open('parsed-story-of-the-stone.txt', 'w', encoding='utf-8') as chineseout, \
     open('utf-story-of-the-stone.txt', 'w', encoding='utf-8') as utfout:
    for line in f:
        for c in line:
            if '\u4e00' <= c <= '\u9fff':  # basic Chinese block
                chineseout.write(c + '\n')
                utfout.write(f"U+{ord(c):04X}\n")

        
