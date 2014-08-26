#!/usr/bin/env python3
# -*- coding: utf-8 -*-
## Usage: python3 task_semcor_wsd.py > task/semcor_wsd.txt

## Extract sentences for different senses of lemmas specified below from SemCor
targets = (("capital", "n"), ("interest", "n"), ("motion", "n"), ("plant", "n"), ("space", "n"), ("suit", "n"), ("tank", "n"), ("vessel", "n")) # natural ambiguous nouns from Schuetze (1998)
## "ruling" (2nd sense = verb gerund) and "train" (noun vs. verbs) have been excluded
targets += (("bank", "n"), ("hand", "n"), ("room", "n")) # some other interesting ambiguous nouns
targets += (("find", "v"), ("grasp", "v"), ("open", "v"), ("run", "v")) # try some verbs

import sys
import nltk
from nltk.corpus import semcor, wordnet
from nltk.stem.wordnet import WordNetLemmatizer
wnl = WordNetLemmatizer()

files = [id for id in semcor.fileids() if not id.startswith("brownv")]
# files = "brown1/tagfiles/br-j05.xml" # for testing

sents_sem = semcor.tagged_sents(fileids=files, tag="sem")
sents_pos = semcor.tagged_sents(fileids=files, tag="pos")
# sents_sem = sents_sem[10:13] # for testing
# sents_pos = sents_pos[10:13]

pos2wn = {"N": wordnet.NOUN, "V": wordnet.VERB, "J": wordnet.ADJ, "R": wordnet.ADV}

print("id\tsid\ttarget\tpos\tsense\tgloss\tsentence\thw\tlemma")
s_num = 0
item_num = {t: 0 for t in targets}
for s_sem, s_pos in zip(sents_sem, sents_pos):
    s_num += 1
    
    lemmas = [x.label() for x in s_sem if isinstance(x, nltk.tree.Tree)] # annotated WordNet lemmas (= synset.hw)
    lemmas = [x for x in lemmas if isinstance(x, nltk.corpus.reader.wordnet.Lemma)] # skip entries where sense isn't a Lemma object

    lemma_info = dict()
    for lemma in lemmas:
        hw = lemma.name()
        synset = lemma.synset()
        pos = synset.name().split(".")[1]
        entry = (hw, pos)
        if entry in lemma_info:
            lemma_info[entry].append(synset)
        else:
            lemma_info[entry] = [synset]
    # print(lemma_info)

    found_targets = [x for x in targets if x in lemma_info]
    if len(found_targets) > 0:

        tokens = [(y, x.label()) for x in s_pos for y in x.leaves()] # distribute POS tags to individual tokens
        pos = [x[1][0] if x[1] is not None else "S" for x in tokens]
        words = [x[0] for x in tokens]
        hws = [wnl.lemmatize(w, pos=pos2wn[p]).lower() if p in pos2wn else w.lower() for w, p in zip(words, pos)]
        hw_pos = ["%s_%s" % x for x in zip(hws, pos)]
        sent_word = " ".join(words)
        sent_hw = " ".join([hw for hw, p in zip(hws, pos) if p != "S"])
        sent_hw_pos = " ".join([x for x in hw_pos if not x.endswith("S")])
        # print(pos)
        # print(hws)

        items = []
        for target in found_targets:
            synsets = set(lemma_info[target])
            if len(synsets) > 1:
                print(target, "is ambiguous in sentence (skipped)", file=sys.stderr, flush=True)
            else:
                sense = synsets.pop()
                item_num[target] += 1
                id = "%s.%s_%d" % (target[0], target[1], item_num[target])
                print("\t".join([id, "%04d" % s_num, target[0], target[1], sense.name(), sense.definition(), sent_word, sent_hw, sent_hw_pos]), flush=True)
