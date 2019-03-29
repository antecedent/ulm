import re
import math
import sys
from collections import defaultdict, Counter
from random import sample
from itertools import combinations

SEPARATOR_PATTERN = re.compile(r'[\W\d_]+', re.UNICODE)
FREQUENCY_CUTOFF = 3
MAX_FREQUENCY_DROP = 10
NUM_COHORTS = 35
COHORT_DEGREES_OF_FREEDOM = 16
MAX_COHORT_MENDINGS = 5
COHORT_CONSENSUS = 7
COHORT_SAMPLE_SIZE = 10000

def make_prefix_tree(words):
    tree = {'_count': 0}
    for word in words:
        node = tree
        node['_count'] += words[word]
        for letter in (word + '#'):
            if letter not in node:
                node[letter] = {'_count': 0}
            node = node[letter]
            node['_count'] += words[word]
    return tree

def valid_child_for_signature(node, child):
    return MAX_FREQUENCY_DROP * node[child]['_count'] > node['_count']        and node[child]['_count'] > FREQUENCY_CUTOFF        and (child == '#' or '#' in node[child]['_signature'])

def collect_signatures(node, nodes_by_signature):
    children = [k for k in node if k[0] != '_']
    for letter in children:
        if letter[0] == '_':
            continue
        collect_signatures(node[letter], nodes_by_signature)
    top_children_by_frequency = sorted(children, key=lambda k: node[k]['_count'], reverse=True)
    top = [k for k in top_children_by_frequency if valid_child_for_signature(node, k)]
    node['_signature'] = ' '.join(sorted('(' + t + ' ' + node[t]['_signature'] + ')' for t in top[:2]))
    if len(top[:2]) == 2:
        nodes_by_signature[node['_signature']].append(node)

def collect_legal_affixes_from(node, legal_affixes, prefix=''):
    if '#' in prefix:
        legal_affixes[prefix[:-1]] += 1
    for letter in node:
        if letter[0] == '_':
            continue
        collect_legal_affixes_from(node[letter], legal_affixes, prefix + letter)

def collect_legal_affixes(signatures, tree, prefix=''):
    legal_affixes = Counter()
    for s in signatures:
        for node in nodes_by_signature[s]:
            collect_legal_affixes_from(node, legal_affixes)
    return set(s for (s, f) in legal_affixes.most_common(int(len(legal_affixes) / math.log(len(legal_affixes)))))

def generate_cohort_index(words):
    index = defaultdict(set)
    for word in words:
        if len(word) > 3:
            for i in range(len(word) - 3):
                for j in range(i + 1, len(word) - 2):
                    for k in range(j + 1, len(word) - 1):
                        for l in range(k + 1, len(word)):
                            index[word[i] + word[j] + word[k] + word[l]].add(word)
    return index

def get_cohort(word, cohort_index, analyses):
    result = set()
    for i in range(len(word) - 3):
        for j in range(i + 1, len(word) - 2):
            for k in range(j + 1, len(word) - 1):
                for l in range(k + 1, len(word)):
                    result |= cohort_index[word[i] + word[j] + word[k] + word[l]]
    result = [word] + list(sample(result, 20))
    while '#'.join(map(lambda c: analyses[c], result)).count(' ') + sum(map(lambda c: math.log(len(c) - 1, 2), result)) > COHORT_DEGREES_OF_FREEDOM:
        result.pop()
    return result

def analyze(node, A, prefix=''):
    if '#' in prefix:
        prefix = prefix[:-1].strip()
        A[prefix.replace(' ', '')] = prefix
    for letter in node:
        if letter[0] == '_':
            continue
        analyze(node[letter], A, prefix + (' ' if '_cut' in node else '') + letter)

def merge_analyses(analyses, reverse_analyses):
    for word in set(analyses):
        a = analyses[word]
        b = reverse_analyses[''.join(reversed(word))]
        consensus = ''
        i = 0
        j = len(b) - 1
        while i < len(a):
            letter = a[i]
            if b[j] == letter:
                consensus += letter
                i += 1
                j -= 1
            else:
                if letter == ' ':
                    consensus += ' '
                    i += 1
                elif b[j] == ' ':
                    consensus += ' '
                    j -= 1
                else:
                    pass
                    #print('???')
        analyses[word] = consensus

print('Usage: python segment.py CORPUS_FILE')

print('Reading corpus...')
raw_corpus = open(sys.argv[1], encoding='utf-8').read().lower()
raw_reverse_corpus = ''.join(reversed(raw_corpus))

all_words = Counter(re.split(SEPARATOR_PATTERN, raw_corpus))
all_reverse_words = Counter(re.split(SEPARATOR_PATTERN, raw_reverse_corpus))

nodes_by_signature = defaultdict(list)
nodes_by_reverse_signature = defaultdict(list)

analyses = {}
reverse_analyses = {}

print('Making forward prefix tree...')
normal_prefix_tree = make_prefix_tree(all_words)

print('Making reverse prefix tree...')
reverse_prefix_tree = make_prefix_tree(all_reverse_words)


print('Collecting forward signatures...')
collect_signatures(normal_prefix_tree, nodes_by_signature)

print('Collecting reverse signatures...')
collect_signatures(reverse_prefix_tree, nodes_by_reverse_signature)

top_signatures = sorted(nodes_by_signature.keys(), key=lambda k: len(nodes_by_signature[k]), reverse=True)
top_reverse_signatures = sorted(nodes_by_reverse_signature.keys(), key=lambda k: len(nodes_by_reverse_signature[k]), reverse=True)

relevant = top_signatures[:int(len(top_signatures) / math.log(len(top_signatures)))]
relevant_reverse = top_reverse_signatures[:int(len(top_reverse_signatures) / math.log(len(top_reverse_signatures)))]

relevant = [k for k in relevant if len(nodes_by_signature[k]) >= FREQUENCY_CUTOFF]
relevant_reverse = [k for k in relevant_reverse if len(nodes_by_reverse_signature[k]) >= FREQUENCY_CUTOFF]

for key in relevant:
    for node in nodes_by_signature[key]:
        node['_cut'] = True

for key in relevant_reverse:
    for node in nodes_by_reverse_signature[key]:
        node['_cut'] = True

legal_suffixes = collect_legal_affixes(relevant, normal_prefix_tree)
legal_prefixes = set(collect_legal_affixes(relevant_reverse, reverse_prefix_tree))

print('Generating cohort index...')
cohort_index = generate_cohort_index(set(w for w, f in all_words.most_common(COHORT_SAMPLE_SIZE)))

print('Separating suffixes...')
analyze(normal_prefix_tree, analyses)

print('Separating prefixes...')
analyze(reverse_prefix_tree, reverse_analyses)

print('Merging suffixal and prefixal analyses...')
merge_analyses(analyses, reverse_analyses)

word = input('> ').strip()
while word != '':
    if word in analyses:
        try:
            vote = Counter()
            for iteration in range(NUM_COHORTS):
                C = get_cohort(word, cohort_index, analyses)
                for i, cw in enumerate(C):
                    C[i] = '#'.join(map(lambda a: '+'.join(list(a)), analyses[cw].split(' ')))
                C = ';'.join(C)
                best_analysis = None
                best_score = None
                gaps = [i for i, c in enumerate(C) if c == '#']
                for count in range(min(MAX_COHORT_MENDINGS, len(gaps)) + 1):
                    for combo in combinations(gaps, count):
                         for addition in [None] + [i for i, c in enumerate(C) if c == '+']:
                            for gap in combo:
                                analysis = C
                                if addition is not None:
                                    analysis = (analysis[:addition] + '#' + analysis[addition + 1:])
                                analysis = (analysis[:gap] + '+' + analysis[gap + 1:])
                                ok = True
                                for w in analysis.replace('+', '').split(';')[:1]:
                                    if len(set(w.split('#'))) < len(w.split('#')):
                                        ok = False
                                        break
                                if not ok:
                                    continue
                                for some_word in analysis.replace('+', '').split(';')[:1]:
                                    morphs = some_word.split('#')
                                    for j in range(1, len(morphs)):
                                        non_stem = ''.join(morphs[j:])
                                        if non_stem not in legal_suffixes:
                                            if ok:
                                                ok = ''.join(morphs[:j]) in legal_prefixes
                                if not ok:
                                    continue
                                item_list = list(map(lambda x: x.replace('+', ''), re.split('[;#]', analysis)))
                                items = Counter(item_list)
                                score = sum(math.log(len(item_list) / k, 2) for (j, k) in items.most_common()) + sum(len(item) * math.log(27, 2) for item in set(items))
                                if  best_score is None or score < best_score:
                                    best_analysis = analysis
                                    best_score = score
                if best_analysis is not None:
                    items = best_analysis.split(';')
                    the_word = [item.replace('+', '').replace('#', ' ') for item in items if item.replace('+', '').replace('#', '') == word][0]
                    vote[the_word] += 1
            print(' ', ', '.join(a for a, v in vote.most_common() if v >= COHORT_CONSENSUS or v == vote.most_common(1)[0][0]))
        except ValueError:
            print(' ', analyses[word])
    elif word[-1] == '!' and word[:-1] in analyses:
        print(' ', analyses[word[:-1]])
    else:
        print(' ', word, '(not in corpus)')
    word = input('> ').strip()
