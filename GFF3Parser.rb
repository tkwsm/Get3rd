#!/usr/bin/ruby

require 'rubygems'
require 'bio'
include Bio

#############################################

module GFF3Parser

  def create_gene_hash
    "bar"
  end
  def GFF3Parser.create_gene_hash

    h = {}
    type_array = %w( gene )
    type_array.each do |type|
      h[type] = []
    end

    return h

  end # create_gene_hash

  def GFF3Parser.create_transcript_hash
  
    h = {}
    type_array = %w( CDS exon intron start_codon stop_codon mRNA transcription_end_site transcription_start_site )
    type_array.each do |type|
      h[type] = []
    end

    return h

  end # create_transcript_hash


  def GFF3Parser.sort_features( features_array )
  
    if features_array == []
    elsif features_array[0][6] == "+"
      features_array.sort!{|x, y| x[3].to_i <=> y[3].to_i }
    elsif features_array[0][6] == "-"
      features_array.sort!{|x, y| y[3].to_i <=> x[3].to_i }
    else
      p "something wrong #{features_array[0]} 2"
      p "something wrong #{features_array[0][6]} 3"
      exit
    end
    return features_array

  end # sort_features( features_array )

  def GFF3Parser.parse_gff3( line, ftype, h, counter_number )

    dataline = []
    seqid  = "" ;  
    source = ""; 
    sstart = 0; 
    send   = 0; 
    score = ".";
    strand = "+"; 
    phase  = "."; 
    attributes = "";
    gid    = ""; 
    tid = "";
    type_array = %w( CDS exon intron start_codon stop_codon transcription_end_site transcription_start_site )
  
    a = line.chomp.split("\t")
    seqid  =  a[0]
    source = a[1]
    ftype   = ftype
    sstart = a[3]
    send   = a[4]
    score  = a[5]
    strand = a[6]
    phase  = a[7]
    attributes = a[8]
    dataline = [seqid, source, ftype, sstart, send, score, strand, phase ]
  
    if ftype == "gene" and attributes =~ /ID/
  
      gid = attributes.slice(/ID=([^\;]+)/, 1)
  
      h[gid] = GFF3Parser.create_gene_hash if h[gid] == nil
      h[gid][ ftype ] = dataline + ["ID=#{gid}\;Name=#{gid}"]
  
    elsif ftype == "mRNA" and attributes =~ /Parent/
  
      gid = attributes.slice(/Parent=([^\;]+)/, 1)
      gid = gid
      tid = attributes.slice(/ID=([^\;]+)/, 1)
      tid = tid
  
      h[gid] = GFF3Parser.create_gene_hash if h[gid] == nil
      h[gid][tid] = GFF3Parser.create_transcript_hash if h[gid][tid] == nil
      h[gid][tid][ ftype ] << dataline + ["ID=#{tid}\;Parent=#{gid}"]

    elsif ftype == "mRNA"
  
      gid = gid
      tid = attributes.slice(/ID=([^\;]+)/, 1)
      gid = tid
      tid = tid
  
      h[gid] = GFF3Parser.create_gene_hash if h[gid] == nil
      h[gid][tid] = GFF3Parser.create_transcript_hash if h[gid][tid] == nil
      h[gid][tid][ ftype ] << dataline + ["ID=#{tid}\;Parent=#{gid}"]
  
    elsif type_array.include?( ftype ) and attributes =~ /Parent/
  
      tid = attributes.slice(/Parent=([^\;]+)/, 1)
      gid = tid.slice(/=?(\S+).t\d/, 1)
      gid = gid
      tid =  tid
  
      h[gid] = GFF3Parser.create_gene_hash if h[gid] == nil
      h[gid][tid] = GFF3Parser.create_transcript_hash if h[gid][tid] == nil
      h[gid][tid][ ftype ] << dataline + ["ID=#{tid}\;Parent=#{gid}"]
  
    else
      p "something bad #{line} 4"
      exit
  
    end
  
    return h 
  
  end # parse_gff3( line, ftype, h, counter_number )

end # module GFF3Parser

###########################################################

codon_h = CodonTable[1].table
codon_rh = {}
codon_h.keys.each do |key|
  v = codon_h[ key ]
  codon_rh[ v ] = [] if codon_rh[ v ] == nil
  codon_rh[ v ] << key
end
# codon_rh.each_key{ |k| print "#{k}\t#{h[k].size}\n" }

# c = { S=>6, L=>6, P=>4, R=>6, T=>4, V=>4, A=>4, G=>4 }

###########################################################
if $0 == __FILE__
###########################################################

  if ARGV.size == 0
  end

  input_file  = ARGV.shift

  GFF3Parser.main_parser( input_file )

###########################################################
end  # __FILE__
###########################################################
#
