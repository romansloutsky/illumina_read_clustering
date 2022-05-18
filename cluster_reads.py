from __future__ import division
import os
from sys import stderr,argv
from multiprocessing import cpu_count
from collections import namedtuple,defaultdict,Counter
from cStringIO import StringIO
import itertools
import HTSeq as htseq
from biomolcore.io import fasta,write_to_file_or_string
from projshared.controllers import EMBOSSAlnController as EMBOSScontrollerBase
from cliceo import CommandLineCaller,PoolManager
from cliceo.workerpool import LabeledObject


#===============================================================================
#         Utils
#===============================================================================

ILLUMLEFT = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC'
ILLUMRIGHT = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG'

REFCOMPDICT = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
def revcomp(rawstr):
  return ''.join([REFCOMPDICT[c] for c in rawstr[::-1]])

def identfrac(zipped):
  num_identities = 0
  for c1,c2 in zipped:
    if c1 == c2:
      num_identities += 1
  return num_identities/len(zipped)

#------------------------------------------------------------------------------ 
# Debugging Utils

def show_overlap(overlap):
  localoverlap = [p for p in overlap]
  while localoverlap:
    if len(localoverlap) > 80:
      front = []
      for i in xrange(80):
        front.append(localoverlap.pop(0))
      print ''.join([c1 for c1,c2 in front])
      print ''.join(['|' if c1==c2 else ' ' for c1,c2 in front])
      print ''.join([c2 for c1,c2 in front])+'\n'
    else:
      print ''.join([c1 for c1,c2 in localoverlap])
      print ''.join(['|' if c1==c2 else ' ' for c1,c2 in localoverlap])
      print ''.join([c2 for c1,c2 in localoverlap])
      localoverlap = False

def determine_and_show_overlap(Sread,Fread,Rread):
  '''
  *** This function is incomplete! ***
  '''
  assert len(Fread) == len(Rread)
  assert len(Sread) <= len(Fread)+len(Rread)
  needlealn = align_unlabeled_str(Fread,revcomp(Rread))
  charpairs = zip(needlealn.seqs[0],needlealn.seqs[1])
  frontgap = []
  while True:
    if '-' in charpairs[0]:
      frontgap.append(charpairs.pop(0))
    else:
      break
  reargap = []
  while True:
    if '-' in charpairs[-1]:
      reargap.append(charpairs.pop())
    else:
      break
  show_overlap(charpairs)


#===============================================================================
#         Fasta and FastQ reads IO
#===============================================================================

#------------------------------------------------------------------------------ 
#         FastQ

def loadFastQ(fpath,raw_iterator=False,**kwargs):
  '''
  Parameters:
    fpath is the path to a FastQ file

    raw_iterator=False (default) => returned items are SequenceWithQualities objects
    raw_iterator=True => returned items are tuples of (seq,name,qualsymbols,qualscale)

    qual_scale has 3 allowed values: 'phred' (default), 'solexa', and 'solexa-old'
  '''
  FQreader = htseq.FastqReader(fpath,raw_iterator=raw_iterator,**kwargs)
  return [r for r in FQreader]


def loadPairedFastQ(fpathF,fpathR,raw_iterator=False,**kwargs):
  '''
  Parameters:
    fpathF is the path to the FastQ file containing F reads
    fpathR is the path to the FastQ file containing R reads

    See loadFastQ() for the rest
  '''
  FQreaderF = htseq.FastqReader(fpathF,raw_iterator=raw_iterator,**kwargs)
  Freads = [r for r in FQreaderF]
  FQreaderR = htseq.FastqReader(fpathR,raw_iterator=raw_iterator,**kwargs)
  Rreads = [r for r in FQreaderR]
  for i,(f,r) in enumerate(zip(Freads,Rreads)):
    if isinstance(f,tuple):
      fname = f[1]
      rname = r[1]
    else:
      fname = f.name
      rname = r.name
    try:
      assert fname.split()[0] == rname.split()[0]
    except AssertionError:
      raise ValueError(
        'In read pair %d, Fread name %s != Rread name %s' % (i,fname.split()[0],
                                                               rname.split()[0]))
    return zip(Freads,Rreads)

def _get_read_parts(readF,readR):
  if isinstance(readF,tuple):
    label = readF[1].split()[0]
    seqF = readF[0]
    seqR = readR[0]
  else:
    label = readF.name.split()[0]
    seqF = readF.seq
    seqR = readR.seq
  return label,seqF,seqR

def _write_read_pair(wh,header,seqF,seqR):
  wh.write('>'+header % 'F'+'\n')
  wh.write(seqF+'\n')
  wh.write('>'+header % 'R'+'\n')
  wh.write(seqR+'\n')

def writePairedFasta(fpath,reads,number=True):
  with open(fpath,'w') as wh:
    if number:
      for i,(f,r) in enumerate(reads):
        label,seqF,seqR = _get_read_parts(f,r)
        label = str(i+1)+'%s'
        _write_read_pair(wh,label,seqF,seqR)
    else:
      for f,r in reads:
        label,seqF,seqR = _get_read_parts(f,r)
        full_label = label.split()[0]+' %s'
        _write_read_pair(wh,full_label,seqF,seqR)

#------------------------------------------------------------------------------ 
#         Fasta

@write_to_file_or_string('target')
def fastawrite_base(label,seq,target=None):
  target.write('>'+label+'\n')
  target.write(seq+'\n')

#------------------------------------------------------------------------------ 
# Single reads

Read = namedtuple('Read',['label','dir','seq'])

def Reads_fromFastQ(fastqreads,label=None):
  origlabel,seqF,seqR = _get_read_parts(*fastqreads)
  label = origlabel if label is None else label
  if not isinstance(label,basestring):
    label = str(label)
  return Read(label,'F',seqF),Read(label,'R',seqF)

def load_reads(src):
  intermediate = fasta.load(src)
  return [Read(r.id[:-1],r.id[-1],str(r.seq)) for r in intermediate]

@write_to_file_or_string('target')
def write_a_read(read,target=None):
  fastawrite_base(read.label+read.dir,read.seq,target=target)

#------------------------------------------------------------------------------ 
# Read pairs

ReadPair = namedtuple('ReadPair',['label','seqF','seqR'])

def load_paired_reads(src):
  intermediate = load_reads(src)
  paired_reads = []
  for i in xrange(0,len(intermediate),2):
    f = intermediate[i]
    r = intermediate[i+1]
    if not f.label == r.label:
      ValueError("Paired reads must be written sequentially with matching "\
                 "labels. Use load_reads() to load unpaired reads.")
    paired_reads.append(ReadPair(f.label,f.seq,r.seq))
  return paired_reads

@write_to_file_or_string('target')
def write_read_pair(readpair,target=None,userevcomp=True):
  target.write('>'+readpair.label+'F\n')
  target.write(readpair.seqF+'\n')
  target.write('>'+readpair.label+'R\n')
  if userevcomp:
    target.write(revcomp(readpair.seqR)+'\n')
  else:
    target.write(readpair.seqR+'\n')
 
@write_to_file_or_string('target')
def write_reads(reads,target=None,direction=None):
  '''
  Parameters:
    :reads is the collection of reads or read pairs to be written
     
    :target is the path to the reads file to be written. If set to None,
      a string is returned instead
   
    :direction determines which subset of reads to write when a collection of
      read pairs is supplied. If set to None, all reads are written. Ignored if
      reads, not read pairs, are supplied.
  '''
  if not all(isinstance(r,Read) for r in reads) and\
                                not all(isinstance(r,ReadPair) for r in reads):
    raise ValueError('reads must be supplied as a uniform collection of Read '\
                     ' or ReadPair instances')
  elif all(isinstance(r,ReadPair) for r in reads):
    if direction == 'F':
      reads_to_write = [Read(r.label,'F',r.seqF) for r in reads]
    elif direction == 'R':
      reads_to_write = [Read(r.label,'R',r.seqR) for r in reads]
    elif direction is None:
      reads_to_write = []
      for read in reads:
        f = Read(read.label,'F',read.seqF)
        r = Read(read.label,'R',read.seqR)
        reads_to_write.append(f)
        reads_to_write.append(r)
    else:
      raise ValueError('direction must be F, R, or None')
  else:
    reads_to_write = reads
   
  for read in reads_to_write:
    write_a_read(read,target)


#===============================================================================
#         Controllers
#===============================================================================


def _passthrough_seqwriter(seqstr):
  '''
  Must be a module-level function because lambda functions and class methods
  can't be pickled, meaning they can't be passed back by worker pool workers.
  
  Probably should rethink this aspect of worker pool design.
  '''
  return seqstr

class EMBOSScontroller(EMBOSScontrollerBase):
  '''
  Since worker pools accept a single argument for each item to be mapped, this
  class exists to handle pairs of strings to align passed as single argument.
  '''
  
  def __init__(self,seqs,EMBOSSprogram,seqwriter=None,write_output_to=None,
                    output_format='srspair',**kwargs):
    seq1,seq2 = seqs
    seqwriter = _passthrough_seqwriter if seqwriter is None else seqwriter
    EMBOSScontrollerBase.__init__(self,EMBOSSprogram,seq1,seq2,seqwriter,
                                  write_output_to=write_output_to,
                                  output_format=output_format,**kwargs)
  
  def call(self):
    self.result = EMBOSScontrollerBase.call(self)


class LinsiController(CommandLineCaller):
  def __init__(self,seqs,seqwriter=fasta.write,capture_stdout=True,
                    silence_stderr=True,**kwargs):
    self.seqs = seqs
    self.seqwriter = seqwriter
    CommandLineCaller.__init__(self,'linsi',capture_stdout=capture_stdout,
                                    silence_stderr=silence_stderr,**kwargs)
  
  def call(self):
    seqs_fname = self.cliCM.write_to_tempfile(self.seqwriter(self.seqs),
                                              dirpath='.')
    self.callstr += ' %s' % seqs_fname
    CommandLineCaller.call(self)
    if self.captured_stdout is not None:
      result = fasta.load(self.captured_stdout)
    else:
      result = None
    return result


class FLAShController(CommandLineCaller):
  def __init__(self,readsfile1,readsfile2,min_overlap=10,max_mismatch_density=0.25,
               max_overlap=251,threads=1):
    callstr = ('flash2 --allow-outies --no-discard true --min-overlap=%d '\
               '--max-overlap=%d --max-mismatch-density=%f --threads=%d '\
               '%s %s') % (min_overlap,max_overlap,max_mismatch_density,threads,
                           os.path.abspath(readsfile1),os.path.abspath(readsfile2))
    CommandLineCaller.__init__(self,callstr,err_to_out=True,capture_stdout=True,
                               in_tmpdir=True,tmpdir_loc='.')
  
  def call(self):
    CommandLineCaller.call(self)
    try:
      self.stitched = loadFastQ('out.extendedFrags.fastq')
      self.unmet = loadPairedFastQ('out.notCombined_1.fastq',
                                   'out.notCombined_2.fastq')
    except IOError:
      raise ValueError('\nSomething went wrong with the FLASh run. '\
                'Here is the FLASh output:\n'+'-'*80+'\n'+self.captured_stdout)


#===============================================================================
#         Alignment
#===============================================================================

#------------------------------------------------------------------------------ 
# Pairwise alignment with EMBOSS

def align_unlabeled_str(s1,s2,EMBOSStool='needle',**kwargs):
  return EMBOSScontrollerBase(EMBOSStool,'>1\n%s' % s1,'>2\n%s' % s2,
                              _passthrough_seqwriter,**kwargs)()

def multiproc_align_unlabeled_str(sequence_of_pairs,numproc,EMBOSStool='needle',
                                  **kwargs):
  for caller in PoolManager(EMBOSScontroller,sequence_of_pairs,numproc=numproc,
                            EMBOSSprogram=EMBOSStool):
    if isinstance(caller,EMBOSScontroller):
      yield caller.result
    else:
      label,caller = caller
      yield label,caller.result

def align_reads(read1,read2,EMBOSStool='needle',**kwargs):
  return EMBOSScontrollerBase(EMBOSStool,read1,read2,write_a_read,**kwargs)()


#------------------------------------------------------------------------------ 
# Multiple sequence alignment with MAFFT L-INS-i (linsi)

@write_to_file_or_string('target')
def write_tupleseqs(tupleseqs,target=None):
  for label,seq in tupleseqs:
    if not isinstance(label,basestring):
      label = str(label)
    fastawrite_base(label,seq,target)

def linsialign(seqs,seqwriter=fasta.write,**kwargs):
  return LinsiController(seqs,seqwriter=seqwriter,**kwargs)()


#===============================================================================
#         Joining overlapping read pairs
#===============================================================================

#-------------------------------------------------------------------------------
# Use the dedicated tool FLASh to stitch reads. It does 95% of the job.

def stitch_with_flash(readsfile1,readsfile2,min_overlap=10,max_mismatch_density=0.25,
                      max_overlap=251,threads=1):
  caller = FLAShController(readsfile1,readsfile2,min_overlap,max_mismatch_density,
                           max_overlap,threads)
  caller()
  return caller.stitched,caller.unmet,caller.captured_stdout

#-------------------------------------------------------------------------------
# Previous, homegrown solution to read stitching, adds a few more to FLASh results

def find_overlap_wo_aln(Fread,Rread,min_ident_frac=0.8,overlap_gap_coff=0.05):
  if Rread.seq[:5] in Fread.seq:
    from_loc_in_F_to_end_of_F = Fread.seq[Fread.seq.find(Rread.seq[:5]):]
    that_len_from_start_of_R = Rread.seq[:len(from_loc_in_F_to_end_of_F)]
    zipped = zip(from_loc_in_F_to_end_of_F,that_len_from_start_of_R)
    if from_loc_in_F_to_end_of_F == that_len_from_start_of_R:
      return zipped
    else:
      lenoverlap = len(zipped)
      IDfrac = identfrac(zipped)
      if lenoverlap >= 8 and IDfrac >= 0.75:
        return zipped
      elif lenoverlap in {6,7} and IDfrac > 0.8:
        return zipped
      else:
        return
#         needlealn = align_unlabeled_str(from_loc_in_F_to_end_of_F,
#                                         that_len_from_start_of_R)
#         if needlealn.header['Gaps'] <= overlap_gap_coff:
#           return zip(needlealn.seqs[0],needlealn.seqs[1])
#         else:
#           return
#   needlealn = align_reads(Fread,Rread)
#   charpairs = zip(needlealn.seqs[0],needlealn.seqs[1])
#   frontgap = []
#   while True:
#     if '-' in charpairs[0]:
#       frontgap.append(charpairs.pop(0))
#     else:
#       break
#   reargap = []
#   while True:
#     if '-' in charpairs[-1]:
#       reargap.append(charpairs.pop())
#     else:
#       break
#   if len([p for p in charpairs if '-' in p])/len(charpairs) <= overlap_gap_coff:
#     return charpairs
#   else:
#     return

InsertIsOverlap = namedtuple('InsertIsOverlap',['label','left','overlap','right'])
def stitch_read_pair(Fread,Rread,overlap):
  o1seq = ''.join([c1 for c1,_ in overlap if c1 != '-'])
  o2seq = ''.join([c2 for _,c2 in overlap if c2 != '-'])
  try:
    # FLASh calls this, most common class of overlaps "innies"
    assert Rread.seq.find(o1seq) == 0 or Rread.seq.find(o2seq) == 0
    joined = Fread.seq+Rread.seq[len(overlap):]
    return Read(Fread.label,'S',joined)
  except AssertionError:
    try:
      # FLASh calls this class of overlaps "outies"
      assert Fread.seq.find(o1seq) == 0 or Fread.seq.find(o2seq) == 0
      overlapseq = Fread.seq[:len(overlap)]
      right = Fread.seq[len(overlap):]
      pos = Rread.seq.find(o1seq)
      if pos == -1:
        pos = Rread.seq.find(o2seq)
      left = Rread.seq[:pos]
      return InsertIsOverlap(Fread.label,left,overlapseq,right)
    except:
      raise

def does_pair_overlap(Fread,Rread,overlap_gap_coff=0.2):
  overlap = find_overlap_wo_aln(Fread,Rread,overlap_gap_coff=overlap_gap_coff)
  if overlap is None:
    return Fread,Rread
  else:
    return stitch_read_pair(Fread,Rread,overlap)

#-------------------------------------------------------------------------------
# Perform read stitching

def gen_read_pairs(paired_reads_collection,makerevcomp=True):
  for pair in paired_reads_collection:
    Fread = Read(pair.label,'F',pair.seqF)
    if makerevcomp:
      Rread = Read(pair.label,'R',revcomp(pair.seqR))
    else:
      Rread = Read(pair.label,'R',pair.seqR)
    yield Fread,Rread
  
def stitch_overlapping_pairs(Freadsfile,Rreadsfile,flash_threads=1,
                             flash_min_overlap=8,min_ident_frac=0.8,
                             makerevcomp=True,verbose=True,shout_every=1000):
  if verbose:
    print >>stderr,"Attempting to stitch paired end reads"
    print >>stderr,"\tFirst-pass stitching with FLASh2"
  stitched,unmet_by_flash,flash_output = stitch_with_flash(Freadsfile,Rreadsfile,
                                                   min_overlap=flash_min_overlap,
                                           max_mismatch_density=1-min_ident_frac,
                                                           threads=flash_threads)
  allreads = {f.name:(f,r) for f,r in loadPairedFastQ(Freadsfile,Rreadsfile)}
  stitched_by_flash = stitched
  stitched = []
  stitched_read_pairs = []
  for i,sr in enumerate(stitched_by_flash):
    stitched.append(Read(str(i+1),'S',sr.seq))
    f,r = allreads[sr.name]
    stitched_read_pairs.append(ReadPair(str(i+1),f.seq,r.seq))
  num_stitched_by_flash = len(stitched)
  num_total = len(allreads)
  unmet_by_flash = [ReadPair(str(i+num_stitched_by_flash+1),f.seq,r.seq)
                    for i,(f,r) in enumerate(unmet_by_flash)]
  if unmet_by_flash and verbose:
    print >>stderr,"\tFLASh2 stitched %d of %d reads" % (num_stitched_by_flash,
                                                                     num_total)
    print >>stderr,"\tWorking on remaining %d reads\n\t" % len(unmet_by_flash),
  unmet = []
  outies = []
  bad = []
  WTF = []
  for i,p in enumerate(gen_read_pairs(unmet_by_flash,makerevcomp=makerevcomp)):
    if verbose and (i+1) % shout_every == 0:
      print >>stderr,i+1,
    try:
      result = does_pair_overlap(p[0],p[1],0.01)
      if isinstance(result,InsertIsOverlap):
        outies.append(result)
      elif isinstance(result,Read):
        stitched.append(result)
        f,r = p
        stitched_read_pairs.append(ReadPair(f.label,f.seq,r.seq))
      else:
        unmet.append(p)
    except AssertionError:
      bad.append(p)
    except:
      WTF.append(p)
  if verbose:
    additionally_stitched = len(stitched) - num_stitched_by_flash
    print >>stderr,'\n\tStitched additional %d reads for a total of %d out of'\
                   ' %d' % (additionally_stitched,len(stitched),num_total)
    print >>stderr,'\t%d outies' % len(outies)
    if bad or WTF:
      print >>stderr,"\t%d bad and %d WTF reads" % (len(bad),len(WTF))
    print >>stderr,'\t%d reads remain unmet' % len(unmet)
  stitched.extend([Read(read.label,'S',read.overlap) for read in outies])
  outie_read_labels = {r.label for r in outies}
  stitched_read_pairs.extend(r for l,r in allreads.items()
                                                     if l in outie_read_labels)
  bad = [ReadPair(F.label,F.seq,R.seq) for F,R in bad]
  WTF = [ReadPair(F.label,F.seq,R.seq) for F,R in WTF]
  unmet = [ReadPair(F.label,F.seq,R.seq) for F,R in unmet]
  return stitched,stitched_read_pairs,unmet,bad,WTF,flash_output


#===============================================================================
#         Clustering
#===============================================================================

#------------------------------------------------------------------------------ 
# Shortcut: clustering sequences of same length without aligning, then merging
# those clusters by aligning their consensus sequences

def cluster_same_len_reads(same_len_reads,ident_cutoff=0.9,verbose=True):
  '''
  Shortcut: Initially clustering sequences of same length without performing
            pairwise alignment.
  '''
  assert len({len(r.seq) for r in same_len_reads}) == 1
  unfinished = [r for r in same_len_reads]
  rawclusters = []
  total_to_do = len(unfinished)
  if verbose and total_to_do >= 50:
    print >>stderr,"Initial clustering %d reads of length "\
                   "%d:" % (total_to_do,len(same_len_reads[0].seq)),
  while unfinished:
    remain_to_do = len(unfinished)
    unmatched = []
    probe = unfinished.pop(0)
    rawclusters.append([probe])
    for r in unfinished:
      if probe.seq == r.seq or identfrac(zip(probe.seq,r.seq)) >= ident_cutoff:
        rawclusters[-1].append(r)
      else:
        unmatched.append(r)
    unfinished = unmatched
    if verbose and total_to_do >= 50 and len(rawclusters[-1]) > 1:
      print >>stderr, len(rawclusters[-1]),
  uniques = []
  clusters = []
  for c in rawclusters:
    if len(c) > 1:
      clusters.append(c)
    else:
      uniques.extend(c)
  if verbose and total_to_do >= 50:
    print >>stderr,"\t%d singletons remain" % len(uniques)
  return clusters,uniques

Cluster = namedtuple('Cluster',['consensus','size','members'])

def consensus_of_single_len(cluster,freq_cutoff=0.95):
  '''
  Find consensus of cluster of sequencens of same length. Since cluster was
  assembled on the basis of pairwise identity without performing alignments,
  no aligning is necessary here either
  '''
  cluster_size = len(cluster)
  if cluster_size > 2:
    real_cutoff = min(freq_cutoff,1/cluster_size)
    consensus_seq = []
    for i in xrange(len(cluster[0].seq)):
      char_counts = Counter(r.seq[i] for r in cluster).most_common()
      char,count = char_counts[0]
      assert count/cluster_size >= real_cutoff
      consensus_seq.append(char)
    return Cluster(''.join(consensus_seq),cluster_size,[r.label for r in cluster])
  else:
    return Cluster(cluster[0].seq,cluster_size,[r.label for r in cluster])

def gen_cluster_pairs_for_merging(probeclust,clusters):
  for i,clust in enumerate(clusters):
    yield LabeledObject(i,('>probecluster\n%s' % probeclust.consensus,
                           '>cluster\n%s' % clust.consensus))

def merge_single_len_clusters(clusters,ident_cutoff,numalnproc,verbose=True):
  '''
  Merge single-read-length clusters when lengths and consensus sequences are
  nearly identical
  '''
  superclusters = []
  while clusters:
    if verbose:
      print >>stderr,len(clusters),
    unmatched = []
    probecluster = clusters.pop(0)
    superclusters.append([probecluster])
    need_to_align = []
    for clust in clusters:
      if abs(len(probecluster.consensus)-len(clust.consensus)) <= 2:
        need_to_align.append(clust)
#         aln = align_unlabeled_str(probecluster.consensus,clust.consensus)
#         if aln.header['Gaps'].abs <= 2 and aln.header['Identity'].frac >= ident_cutoff:
#           superclusters[-1].append(clust)
#         else:
#           unmatched.append(clust)
      else:
          unmatched.append(clust)
    gen = gen_cluster_pairs_for_merging(probecluster,need_to_align)
    to_supercluster = []
    to_unmatched = []
    for clustnum,aln in multiproc_align_unlabeled_str(gen,numalnproc):
      if aln.header['Gaps'].abs <= 2 and aln.header['Identity'].frac >= ident_cutoff:
        to_supercluster.append(need_to_align[clustnum])
      else:
        to_unmatched.append(need_to_align[clustnum])
    superclusters[-1].extend(to_supercluster)
    unmatched.extend(to_unmatched)
    clusters = unmatched
  return superclusters

#------------------------------------------------------------------------------ 
# Attempt to place singleton reads into clusters

def gen_unique_cluster_pairs_to_aln(uniques,reduced_clusters,minsize,verbose=True,
                                    shout_every=100):
  for i,r in enumerate(uniques):
    if verbose and (i+1) % shout_every == 0:
      print >>stderr,i+1,
    for j,rc in enumerate(reduced_clusters):
      if sum(c.size for c in rc) >= minsize:
        for k,c in enumerate(rc):
          if abs(len(r.seq)-len(c.consensus)) <= 2:
            yield LabeledObject((('unique_read',i),('supercluster',j),
                                 ('cluster',k)),
                                ('>read\n%s' % r.seq,
                                 '>cluster_consensus\n%s' % c.consensus))
  if uniques and verbose and i > shout_every:
    print >>stderr,'\n',

def find_cluster_homes_for_uniques(reduced_clusters,uniques,numalnproc,minsize,
                                   ident_cutoff=0.9,verbose=True,shout_every=100):
  if verbose:
    print >>stderr,"Finding cluster homes for %d previously unclustered" % len(uniques)
  gen = gen_unique_cluster_pairs_to_aln(uniques,reduced_clusters,minsize,
                                        verbose=verbose,shout_every=shout_every)
  matches = defaultdict(list)
  for label,aln in multiproc_align_unlabeled_str(gen,numalnproc):
    label = dict(label)
    ureadnum = label['unique_read']
    superclusternum = label['supercluster']
    if superclusternum not in matches[ureadnum]:
      if aln.header['Gaps'].abs <= 2 and\
                                   aln.header['Identity'].frac >= ident_cutoff:
        matches[ureadnum].append(superclusternum)
  still_unique = []
  ambiguously_matched = []
  for readnum,matched_to in matches.items():
    r = uniques[readnum]
    if len(matched_to) > 1:
      ambiguously_matched.append((r,matched_to))
    elif len(matched_to) == 1:
      reduced_clusters[matched_to[0]].append(Cluster(r.seq,1,[r.label]))
    else:
      still_unique.append(r)
  return ambiguously_matched,still_unique

#------------------------------------------------------------------------------ 
# Weighted consensus of cluster consensuses by L-INS-i alignment

# def full_reads_consensus_w_linsi(reads,warn_if_freq_below=0.6):
#   pass

def weighted_cluster_consensus_w_linsi(supercluster,warn_if_freq_below=0.6):
  aln = linsialign(enumerate([clust.consensus for clust in supercluster]),
                   seqwriter=write_tupleseqs)
  alnseqs = [(int(s.id),str(s.seq).upper()) for s in aln]
  total = sum([clust.size for clust in supercluster])
  real_warn_min_freq = min(warn_if_freq_below,1/total)
  pos_char_freqs = []
  for i in xrange(len(alnseqs[0][1])):
    charcounts = Counter()
    for clustNum,seq in alnseqs:
      char = seq[i]
      charcounts[char] += supercluster[clustNum].size
    pos_char_freqs.append([(char,count/total)
                           for char,count in charcounts.most_common()])
  consensus_seq = []
  for i,cf in enumerate(pos_char_freqs):
    char,freq = cf[0]
    if char != '-':
      consensus_seq.append(char)
      if freq < real_warn_min_freq:
        print >>stderr,'WARNING: at position %d most frequent base is only %0.3f' % (i,freq)
  return Cluster(''.join(consensus_seq),total,[l for c in supercluster
                                                  for l in c.members])

def make_supercluster_consensuses(superclusters):
  sorted_by_size = sorted(superclusters,
                          key=lambda x: sum([c.size for c in x]),reverse=True)
  consensusified = []
  for sc in sorted_by_size:
    if len(sc) > 1:
      consensusified.append(weighted_cluster_consensus_w_linsi(sc))
    else:
      consensusified.append(sc[0])
  return consensusified

#------------------------------------------------------------------------------ 
#  Bundled task: finalize initial single-read-length clusters using
#                pairwise alignment by:
#    a) merging highly similar clusters into superclusters
#    b) attempt to place unique sequences (optional)
#    c) determining consensus sequences for superclusters
#------------------------------------------------------------------------------

def reduce_initial_clusters(clusters,uniques,numalnproc,ident_cutoff=0.9,
                            consensus_freq_cutoff=0.95,match_uniques=True,
                            min_clust_size_for_uniq_matching=10,verbose=True):
  '''
  Bundled task: finalize initial single-read-length clusters using
                pairwise alignment by:
    a) merging highly similar clusters into superclusters
    b) attempt to place unique sequences (optional)
    c) determining consensus sequences for superclusters
  '''
  if verbose:
    print >>stderr,"Reducing %d clusters: consensus ..." % len(clusters),
  reduced_clusters = [consensus_of_single_len(c,consensus_freq_cutoff)
                      for c in clusters]
  if verbose:
    print >>stderr,"merging:",
  superclusters = merge_single_len_clusters(reduced_clusters,ident_cutoff,
                                            numalnproc,verbose=verbose)
  if verbose:# and not match_uniques:
    print >>stderr,"finished"
    print >>stderr,"Merged %d clusters into %d superclusters" % (len(clusters),
                                                                 len(superclusters))
  if match_uniques:
    ambig_matched_uniques,uniques = find_cluster_homes_for_uniques(superclusters,
                                                                   uniques,
                                                                   numalnproc,
                                                  min_clust_size_for_uniq_matching,
                                                                   ident_cutoff,
                                                                   verbose=verbose)
    if not ambig_matched_uniques:
      ambig_matched_uniques = None
  else:
    ambig_matched_uniques = None
  if verbose:
    print >>stderr,"Making supercluster consensuses ...",
  reduced_again = make_supercluster_consensuses(superclusters)
  if verbose:
    print >>stderr,"finished"
  return reduced_again,ambig_matched_uniques,uniques

#------------------------------------------------------------------------------ 
# Cluster reads by one end

OneDirClusterSet = namedtuple('OneDirClusterSet',['dir','clusters','ambiguous',
                                                  'uniques'])

def get_one_dir(paired_reads,direction):
  assert direction in ('F','R')
  if direction == 'F':
    return [Read(p.label,'F',p.seqF) for p in paired_reads]
  elif direction == 'R':
    return [Read(p.label,'R',p.seqR) for p in paired_reads]
  else:
    raise ValueError("direction must be one of %s" % str(('F','R')))

def cluster_one_end(paired_reads,direction,numalnproc,ident_cutoff=0.9,
                    consensus_freq_cutoff=0.95,match_uniques=True,
                    min_clust_size_for_uniq_matching=10,verbose=True):
  paired_reads_are_ReadPairs = all(isinstance(r,ReadPair) for r in paired_reads)
  assert paired_reads_are_ReadPairs or\
         all(isinstance(r,Read) for rp in paired_reads for r in rp)
  if paired_reads_are_ReadPairs:
    reads = get_one_dir(paired_reads,direction)
  else:
    reads = [r for rp in paired_reads for r in rp if r.dir == direction]
  if verbose:
    print >>stderr,"Clustering %s reads: initial clustering" % direction
  clusters,uniques = cluster_same_len_reads(reads,ident_cutoff=ident_cutoff,
                                            verbose=verbose)
  if verbose:
    print >>stderr,"Clustering %s reads: reducing clusters" % direction
  return OneDirClusterSet(direction,*reduce_initial_clusters(clusters,uniques,
                                                             numalnproc,
                                                             ident_cutoff,
                                                         consensus_freq_cutoff,
                                                             match_uniques,
                                              min_clust_size_for_uniq_matching,
                                                             verbose=verbose))

#------------------------------------------------------------------------------ 
# Cluster stitched reads

def cluster_stitched(reads,numalnproc,ident_cutoff=0.9,consensus_freq_cutoff=0.95,
                     match_uniques=True,min_clust_size_for_uniq_matching=10,
                     verbose=True):
  bylength = defaultdict(list)
  for r in reads:
    bylength[len(r.seq)].append(r)
  bylength = sorted(bylength.items(),key=lambda x: len(x[1]),reverse=True)
  clusters = []
  uniques = []
  if verbose:
    print >>stderr,"Clustering %d stitched reads: initial clustering" % len(reads)
  started_on_small_length_batches = False
  for l,same_len_reads in bylength:
    if verbose and l < 50:
      if not started_on_small_length_batches:
        started_on_small_length_batches = True
        print >>stderr, "Finishing up with small read collections:",
      print >>stderr,l,
    c,u = cluster_same_len_reads(same_len_reads,ident_cutoff=ident_cutoff,
                                 verbose=verbose)
    clusters.extend(c)
    uniques.extend(u)
  if verbose:
    if started_on_small_length_batches:
      print >>stderr,'\n',
    print >>stderr,"Clustering %d stitched reads: reducing clusters" % len(reads)
  return OneDirClusterSet('S',*reduce_initial_clusters(clusters,uniques,
                                                       numalnproc,ident_cutoff,
                                                       consensus_freq_cutoff,
                                                       match_uniques,
                                              min_clust_size_for_uniq_matching,
                                                       verbose=verbose))

#------------------------------------------------------------------------------ 
# Cluster reads that don't overlap in the middle by one end,then by the other

PairedClusterSet = namedtuple('PairedClusterSet',['dir','cluster','opp_clusters',
                                                  'opp_ambiguous','opp_uniques'])

BiDirClusterSet = namedtuple('BiDirClusterSet',['firstdir','bidir_clusters',
                                                'opp_clusters','ambiguous',
                                                'opp_ambiguous','bidir_uniques'])

def cluster_unmet(paired_reads,first_direction,numalnproc,ident_cutoff=0.9,
                  consensus_freq_cutoff=0.95,match_uniques=True,
                  min_clust_size_for_uniq_matching=10,verbose=True):
  if verbose:
    print >>stderr,"Clustering %d unmet reads: starting with"\
                   " %s" % (len(paired_reads),first_direction)
  first_dir_clustset = cluster_one_end(paired_reads,first_direction,numalnproc,
                                       ident_cutoff,consensus_freq_cutoff,
                                       match_uniques,min_clust_size_for_uniq_matching,
                                       verbose=verbose)
  opposite_direction = 'F' if first_direction == 'R' else 'R'
  bidir_clusters = []
  if verbose:
    print >>stderr,"Clustering %d unmet reads: clustering %s for each %s"\
                   " cluster" % (len(paired_reads),opposite_direction,
                                 first_direction)
  for clust in first_dir_clustset.clusters:
    if verbose:
      print >>stderr,"Clustering %s for %s cluster of size"\
                     " %d" % (opposite_direction,first_direction,clust.size)
    cluster_members = set(clust.members)
    reads_in_cluster = [r for r in paired_reads if r.label in cluster_members]
    opp_dir_clustset = cluster_one_end(reads_in_cluster,opposite_direction,
                                       numalnproc,ident_cutoff,
                                       consensus_freq_cutoff,match_uniques,
                                       min_clust_size_for_uniq_matching,
                                       verbose=verbose)
    bidir_clusters.append(PairedClusterSet(first_direction,clust,
                                           *opp_dir_clustset[1:]))
  unique_labels = {r.label for r in first_dir_clustset.uniques}
  paired_uniques = [r for r in paired_reads if r.label in unique_labels]
  if verbose:
    print >>stderr,"Clustering %d unmet reads: clustering %s for %d %s "\
                   "unclustered" % (len(paired_reads),opposite_direction,
                                    len(first_dir_clustset.uniques),
                                    first_direction)
  uniques_opp_dir_clustset = cluster_one_end(paired_uniques,opposite_direction,
                                             numalnproc,ident_cutoff,
                                             consensus_freq_cutoff,match_uniques,
                                             min_clust_size_for_uniq_matching,
                                             verbose=verbose)
  bidir_unique_labels = {r.label for r in uniques_opp_dir_clustset.uniques}
  bidir_uniques = [r for r in paired_reads if r.label in bidir_unique_labels]
  return BiDirClusterSet(first_direction,bidir_clusters,
                         uniques_opp_dir_clustset.clusters,
                         first_dir_clustset.ambiguous,
                         uniques_opp_dir_clustset.ambiguous,bidir_uniques)


#===============================================================================
#         Results Output
#===============================================================================

#------------------------------------------------------------------------------ 
# Write results of paired-end read stitching

def write_read_stitching_output(topresultsdir,stitched,stitched_read_pairs,unmet,
                                bad,WTF,flash_output,revcomp_for_R=False):
  if stitched:
    write_reads(stitched,'%s/stitched_reads.fasta' % topresultsdir)
    write_reads(stitched_read_pairs,
                '%s/stitched_read_pairs.fasta' % topresultsdir)
  if unmet:
    write_reads(unmet,'%s/unmet_read_pairs.fasta' % topresultsdir)
  if bad:
    with open('%s/bad.fasta' % topresultsdir,'w') as wh:
      for r in bad:
        write_read_pair(r,wh)#,userevcomp=True)
  if WTF:
    with open('%s/WTF.fasta' % topresultsdir,'w') as wh:
      for r in WTF:
        write_read_pair(r,wh)#,userevcomp=True)
  with open('FLASh_output.txt','w') as wh:
    wh.write(flash_output)

#------------------------------------------------------------------------------ 
# Write clustering results

@write_to_file_or_string('target')
def write_unclustered_reads(reads,target,revcomp_for_R=False):
  stitched = all(isinstance(r,Read) for r in reads)
  paired = all(isinstance(r,ReadPair) for r in reads)
  assert stitched or paired
  if stitched:
    for r in reads:
      write_a_read(r,target)
  else:
    for r in reads:
      write_read_pair(r,target,userevcomp=revcomp_for_R)

def write_unpaired_read_cluster(cluster,rawreads,clusternum,out_of,dirname='.',
                                prefix=''):
  if not os.path.exists(dirname):
    os.mkdir(dirname)
  numdigits = str(len(str(out_of)))
  with open(('%s/%sCluster%0'+numdigits+'dconsensus.fasta') % (dirname,prefix,
                                                             clusternum),
                                                                    'w') as wh:
    wh.write(('>ClustOf%d\n%s\n' % (cluster.size,cluster.consensus)))
  cluster_members = set(cluster.members)
  reads_in_cluster = [r for r in rawreads if r.label in cluster_members]
  with open(('%s/%sCluster%0'+numdigits+'dreads.fasta') % (dirname,prefix,
                                                         clusternum),'w') as wh:
    for r in reads_in_cluster:
      write_a_read(r,wh)

def write_unpaired_read_clusterset(clusterset,rawreads,dirname='.',prefix=''):
  if not os.path.exists(dirname):
    os.mkdir(dirname)
  total_clusters = len(clusterset.clusters)
  for i,clust in enumerate(clusterset.clusters):
    write_unpaired_read_cluster(clust,rawreads,i+1,total_clusters,dirname,prefix)
  if clusterset.ambiguous is not None:
    print >>stderr,'WARNING: Ambiguous reads in %s direction' % clusterset.dir
  if clusterset.uniques:
    write_unclustered_reads(clusterset.uniques,
                            '%s/%sUnclustered.fasta' % (dirname,prefix))

def write_paired_clusterset(clusterset,reads_by_dir,clusternum,out_of,dirname='.'):
  numdigits = str(len(str(out_of)))
  targetdir = ('%s/%sCluster%0'+numdigits+'d') % (dirname,clusterset.dir.lower(),
                                                clusternum)
  if not os.path.exists(targetdir):
    os.makedirs(targetdir)
  write_unpaired_read_cluster(clusterset.cluster,reads_by_dir[clusterset.dir],
                              clusternum,out_of,targetdir,
                              prefix=clusterset.dir.lower())
  oppdir = 'F' if clusterset.dir == 'R' else 'R'
  total_opp_clusters = len(clusterset.opp_clusters)
  for i,oppclust in enumerate(clusterset.opp_clusters):
    write_unpaired_read_cluster(oppclust,reads_by_dir[oppdir],i+1,
                                total_opp_clusters,targetdir,prefix=oppdir+'opp')
  if clusterset.opp_ambiguous is not None:
    print >>stderr,'WARNING: Ambiguous oppdir (%s) reads for cluster '\
                   '%d' % (oppdir,clusternum)
  if clusterset.opp_uniques:
    write_unclustered_reads(clusterset.opp_uniques,
                           '%s/%sUnclustered.fasta' % (targetdir,
                                                       oppdir+'opp'))

def write_bidir_clusterset(clusterset,rawreads,dirname='.'):
  reads_by_dir = {'F':get_one_dir(rawreads,'F'),'R':get_one_dir(rawreads,'R')}
  total_clusters = len(clusterset.bidir_clusters)
  for i,clust in enumerate(clusterset.bidir_clusters):
    write_paired_clusterset(clust,reads_by_dir,i+1,total_clusters,dirname)
  if clusterset.ambiguous is not None:
    print >>stderr,'WARNING: Ambiguous reads in first direction (%s) of bidir'\
                   ' cluster set' % clusterset.firstdir
  oppdir = 'F' if clusterset.firstdir == 'R' else 'R'
  unclustdir = '%s/%sUnclustered' % (dirname,clusterset.firstdir.lower())
  if not os.path.exists(unclustdir):
    os.makedirs(unclustdir)
  total_opp_clusters = len(clusterset.opp_clusters)
  for i,clust in enumerate(clusterset.opp_clusters):
    write_unpaired_read_cluster(clust,reads_by_dir[oppdir],i+1,total_opp_clusters,
                                unclustdir,prefix=oppdir+'opp')
  if clusterset.opp_ambiguous is not None:
    print >>stderr,'WARNING: Ambiguous reads in opposite direction (%s) of bidir'\
                   ' cluster set when clusteriing opp ends of first dir '\
                   'unclustered reads' % oppdir
  if clusterset.bidir_uniques:
    write_unclustered_reads(clusterset.bidir_uniques,
                           '%s/bidirectionallyUnclustered.fasta' % (dirname))


#===============================================================================
#         Top-level Analysis Performers
#===============================================================================


def analyze_one_end_reads(Freadsfile,Rreadsfile,numalnproc,direction='F',
                          ident_cutoff=0.9,consensus_freq_cutoff=0.95,
                          match_uniques=True,min_clust_size_for_uniq_matching=10,
                          topresultsdir='.',verbose=True):
  '''
  If it is a priori known that only one read direction is useful, then only those
  reads need to be clustered
  '''
  
  paired_reads = [Reads_fromFastQ(p,i+1)
                  for i,p in enumerate(loadPairedFastQ(Freadsfile,Rreadsfile))]
  if verbose:
    print >>stderr,"Performing one direction (%s) clustering on %d paired"\
                   " reads" % (direction,len(paired_reads))
  clusterset = cluster_one_end(paired_reads,direction,numalnproc,ident_cutoff,
                               consensus_freq_cutoff,match_uniques,
                               min_clust_size_for_uniq_matching,verbose=verbose)
  targetdir = '%s/%s' % (topresultsdir,'%sreadClustering' % direction)
  if not os.path.exists(targetdir):
    os.makedirs(targetdir)
  oppdir = 'F' if clusterset.dir == 'R' else 'R'
  dirReads = []
  oppReads = {}
  for Fr,Rr in paired_reads:
    tmp = {'F':Fr,'R':Rr}
    dirReads.append(tmp[direction])
    oppReads[tmp[oppdir].label] = tmp[oppdir]
  if verbose:
    print >>stderr,"Writing output ...",
  write_unpaired_read_clusterset(clusterset,dirReads,targetdir,direction.lower())
  numdigits = str(len(str(len(clusterset.clusters))))
  for i,clust in enumerate(clusterset.clusters):
    with open(('%s/%soppCluster%0'+numdigits+'dreads.fasta') % (targetdir,oppdir,i+1),
                                                                    'w') as wh:
      for label in clust.members:
        oppR = oppReads[label]
        wh.write('>%s%s\n%s\n' % (label,oppdir,oppR.seq))
  with open('%s/%soppUnclustered.fasta' % (targetdir,oppdir),'w') as wh:
    for r in clusterset.uniques:
      oppR = oppReads[r.label]
      wh.write('>%s%s\n%s\n' % (r.label,oppdir,oppR.seq))
  if verbose:
    print >>stderr,"Finished"
  

def analyze_read_pairs(Freadsfile,Rreadsfile,numalnproc,first_dir_for_unmet='F',
                       ident_cutoff=0.9,consensus_freq_cutoff=0.95,
                       match_uniques=True,min_clust_size_for_uniq_matching=10,
                       topresultsdir='.',verbose=True):
  if verbose:
    print >>stderr,"Analyzing paired reads from %s, %s" % (Freadsfile,Rreadsfile)
    print >>stderr,"Classifying reads"
  stitched,stitched_read_pairs,unmet,bad,WTF,flash_output =\
       stitch_overlapping_pairs(Freadsfile,Rreadsfile,flash_threads=numalnproc,
                                                               verbose=verbose)
  if verbose:
    print >>stderr,"Writing read stitching results ...",
  write_read_stitching_output(topresultsdir,stitched,stitched_read_pairs,unmet,
                              bad,WTF,flash_output)
  stitched_clusterset = cluster_stitched(stitched,numalnproc,
                                         ident_cutoff=ident_cutoff,
                                   consensus_freq_cutoff=consensus_freq_cutoff,
                                         match_uniques=match_uniques,
              min_clust_size_for_uniq_matching=min_clust_size_for_uniq_matching,
                                         verbose=verbose)
  unmet_clusterset = cluster_unmet(unmet,first_dir_for_unmet,numalnproc,
                                   ident_cutoff=ident_cutoff,
                                   consensus_freq_cutoff=consensus_freq_cutoff,
                                   match_uniques=match_uniques,verbose=verbose)
  if verbose:
    print >>stderr,"Writing output ...",
  for dirname in ['stitchedClustering','unmetClustering']:
    if not os.path.exists('%s/%s' % (topresultsdir,dirname)):
      os.makedirs('%s/%s' % (topresultsdir,dirname))
  write_unpaired_read_clusterset(stitched_clusterset,stitched,
                                 '%s/stitchedClustering' % topresultsdir,
                                 'stitched')
  write_bidir_clusterset(unmet_clusterset,unmet,
                         '%s/unmetClustering' % topresultsdir)
  if verbose:
    print >>stderr,"Finished"

if __name__ == "__main__" and len(argv) >= 3:
  Freadsfile = argv[1]
  Rreadsfile = argv[2]
  numalnproc = int(argv[3])
  targetdir = argv[4]
  first_dir = argv[5] if len(argv) > 5 else 'F'
  verbose = False if 'silent' in argv[6:] else True
  one_dir_only = True if 'onedir' in argv[6:] else False
  match_uniques = False if 'skip_uniques' in argv[6:] else True
  if not os.path.exists(Freadsfile):
    print >>stderr,"%s does not exist, exiting" % Freadsfile
  if not os.path.exists(Rreadsfile):
    print >>stderr,"%s does not exist, exiting" % Rreadsfile
  if not match_uniques:
    print >>stderr,"Will not attempt to find cluster homes for initially unique reads"
  try:
    if not os.path.exists(targetdir):
      os.mkdir(targetdir)
  except:
    print >>stderr,"Could not create directory %s, exiting" % targetdir
  if one_dir_only:
    analyze_one_end_reads(Freadsfile,Rreadsfile,numalnproc,direction=first_dir,
                          match_uniques=match_uniques,topresultsdir=targetdir,
                          verbose=verbose)
  else:
    analyze_read_pairs(Freadsfile,Rreadsfile,numalnproc,
                       first_dir_for_unmet=first_dir,match_uniques=match_uniques,
                       topresultsdir=targetdir,verbose=verbose)

