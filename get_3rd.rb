#!/usr/bin/ruby

require 'rubygems'
require 'bio'
include Bio

class GFFParser

  def initialize( gff3 )
    @gff3genehash  = gff3parser( gff3 )[0]
    @gff3scafhash  = gff3parser( gff3 )[1]
  end

  attr_reader :gff3genehash, :gff3scafhash

  def gff3parser( gff3 )
    gff3genehash = {}   
    gff3scafhash = {}   
######
    a = []
    open(gff3).each do |x|
      a = x.chomp.split("\s")
      fseqid  = a[0] 
      fsource = a[1] 
      ftype   = a[2] 
      fstart  = a[3] 
      fend    = a[4] 
      fscore  = a[5] 
      fstrand = a[6] 
      fattr   = a[7] 
      fid     = fattr.slice(/ID=([^\;]+)/, 1)
      fparent = fattr.slice(/Parent=([^\;]+)/, 1)
      if    ftype == "gene"
        gff3scafhash[ fseqid ] = [] if gff3scafhash[ fseqid ] == nil
        gff3scafhash[ fseqid ] << fid
	if gff3genehash[ fid ] == nil
          gff3genehash[ fid ] = {:cds => [], :strand => "" } 
        end	  
      elsif ftype == "cds"
        gff3genehash[ fid ][:cds]    << [ fstart, fend ]
        gff3genehash[ fid ][:strand] = fstrand
      end
    end
#    gff3genehash = { "gene1" => {:cds => [[7, 12], [22, 30]], 
#	                         :strand => "+"} }
#    gff3scafhash = { "scaf1" => [ "gene1" ] }
    return [ gff3genehash , gff3scafhash ]
  end

end

class FastaParser

  def initialize( fastafile ) 
    @fastahash = fastaparser( fastafile )
  end

  attr_reader :fastahash

  def fastaparser( fastafile )
    fastahash = {}
    FlatFile.new( FastaFormat, open( fastafile ) ).each do |e|
      fastahash[ e.definition ] = e.naseq
    end
    return fastahash
  end

end

class MergeGFFontoFASTA

  def initialize( gff3, fastafile )
    @gff3genehash  = GFFParser.new( gff3 ).gff3genehash
    @gff3scafhash  = GFFParser.new( gff3 ).gff3scafhash
    @fastahash = FastaParser.new( fastafile ).fastahash
  end

  attr_reader :gff3genehash, :gff3scafhash, :fastahash

end

# UTR     "ctataa"
# CDS1    "atg ggc"
# PEP1    "Met Gly"
# INTRON1 "tacgatgct"
# CDS2    "atc ttc taa"
# PEP2    "Ile Phe Stop"
# InterG1 "tagctgatcgatcgatgctagctagctgatcgatg"

if $0 == __FILE__

  gff3      = ARGV.shift
  fastafile = ARGV.shift
  mgf = MergeGFFontoFASTA.new( gff3, fastafile )

  mgf.gff3scafhash.each_key do |scafid|
    mgf.gff3scafhash[ scafid ].each do |gid|
      orf = ""
      mgf.gff3genehash[ gid ][:cds].each do |aregion|
        cdss = aregion[0].to_i
        cdse = aregion[1].to_i
        if    mgf.gff3genehash[ gid ][:strand] == "+" 
          orf += mgf.fastahash[scafid].subseq(cdss, cdse)
        elsif mgf.gff3genehash[ gid ][:strand] == "-" 
          orf += mgf.fastahash[scafid].subseq(cdss, cdse).complement
        end
      end
      if orf.length % 3 != 0
        next STDERR.print "total length of ORF must be multiples of 3!!\n"; 
      end
      Bio::Sequence::NA.new( orf ).window_search(3, 3) do |codon|
        p codon
        p codon.translate
      end
    end
  end

end


