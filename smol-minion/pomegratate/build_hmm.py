from pomegranate import *
import numpy as np
import math
import sys
import pandas
from read_st import read_st_file

def build_hmm(transition_triples, states):
  model = HiddenMarkovModel()

  model.add_states(*states)

  for (r, c, val) in transition_triples:
    source = None
    if r == 0:
      source = model.start
    else:
      source = states[r - 1]

    dest = None
    if c == len(states) + 1:
      dest = model.end
    elif c == 0:
      dest = model.start
    else:
      dest = states[c - 1]

    model.add_transition(source, dest, val)

  print "done triples"
  model.bake()
  print "done bake"
  model.freeze()
  model.freeze_distributions()

  return model

def check(triples, states):
  n_states = len(states)

  rows = [t[0] for t in triples]
  cols = [t[1] for t in triples]

  assert(min(rows) == 0)
  assert(max(cols) == n_states + 1) # -1 because of 1 indexing, +2 because of start, end column

def read_emission_csv(filepath):
  return np.loadtxt(filepath, delimiter=",")

def read_emission_pickle(filepath):
  return pandas.read_pickle(filepath).as_matrix()

def read_st_file(filepath):
  with open(filepath) as f:
    lines = f.readlines()
    tuples = []
    for l in lines:
      if l[0] == '#':
        continue
      words = l.split()
      r = int(words[0])
      c = int(words[1])
      val = float(words[2])
      tuples.append((r, c, val))
    return tuples

def modelIxToStateIx(model, modelIx):
  try:
    return int(model.states[modelIx].name)
  except ValueError:
    return None

def statePermutation(model):
  perm = filter(lambda i: i, map(lambda s: modelIxToStateIx(model, s), range(len(model.states))))
  #iperm = map(lambda i: perm.index(i + 1), range(len(model.states)))

  iperm = sorted(xrange(len(perm)), key=lambda ix: perm[ix])

  return np.array(perm), np.array(iperm)

def build_emissions_hmm(trans_file = "test.st",
                        emissions_file = "test.csv",
                        posterior_file = None,
                        viterbi_file = None,
                        forward_file = None,
                        backward_file = None):
  emissions = read_emission_csv(emissions_file) if emissions_file.endswith(".csv") \
              else read_emission_pickle(emissions_file)

  n_loci, n_states = emissions.shape

  states = [State(DiscreteDistribution({i: pi for (i, pi) in enumerate(state_col)}),
                                       name = str(s))
            for s, state_col in enumerate(emissions.T)]

  triples = read_st_file(trans_file)

  print "read"

  model = build_hmm(triples, states)

  print "built"

  path = [i for i in range(n_loci)]

  perm, iperm = statePermutation(model)
  vit = viterbi_seq(model, path)
  print "done viterbi"
  post = posterior_seq(model, path)[:, iperm]
  print "done post"
  fwd = forward_seq(model, path)
  bwd = backward_seq(model, path)

  if viterbi_file:
    with open(viterbi_file, 'w') as f:
      f.write(filter(lambda x: x!=' ', str(vit)[1:-1]))

  if posterior_file:
    np.savetxt(posterior_file, post, delimiter = ",")

  if forward_file:
    np.savetxt(forward_file, fwd, delimiter = ",")

  if backward_file:
    np.savetxt(backward_file, bwd, delimiter = ",")

  # model.draw() might be cool, need an additional package
  return model, emissions, post, vit, fwd, bwd

def forward_seq(model, path):
  emissions = model.forward(path)
  return math.e ** emissions

def backward_seq(model, path):
  emissions = model.backward(path)
  return math.e ** emissions

def posterior_seq(model, path):
  (t, emissions) = model.forward_backward(path)
  return math.e ** emissions

def viterbi_seq(model, path):
  v = model.viterbi(path)
  return [int(s[1].name) for s in v[1][1:-1]]

if len(sys.argv) > 2:
  build_emissions_hmm(*sys.argv[1:])
