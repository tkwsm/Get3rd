#!/usr/bin/ruby

require 'rubygems'
require 'bio'
include Bio

module ParseAlign

  def huga
    "foo"
  end

  class HandleAlignment

    def initialize( consensus_hash )
      @c = consensus_hash
      @gids  = {}; @c[0].each_with_index{ |gid, i| @gids[ gid ] = i }
      @i2gid = {}; @c[0].each_with_index{ |gid, i| @i2gid[ i ] = gid }
      @a = get_aligned
      @ac = get_all_conserved
      @cp = get_compensation
    end

    attr_reader :c, :a, :ac, :cp, :gids, :i2gid

    def get_compensation
      cp = {}
      # @c: keys is the alley of alignment-position
      alignment_positions = @c.keys.sort
      num_of_species = @c.values[0].size
      alignment_positions.each{|k| cp[k] = Array.new( num_of_species )}
      non_gap_pos = []
      for s in 0..( num_of_species - 1 )
        non_gap_pos = alignment_positions.collect{|k| k if @c[k][s] != "-" and @c[k][s] != "*" }.compact.sort
        non_gap_pos.each_with_index do |k, i|
          cp[k][s] = i
        end
      end
      return cp
    end

    def get_aligned
      aligned = {}
      @c.each_key do |k|
        unless @c[k].include?("-") or @c[k].include?("*") 
          aligned[k] = @c[k]
        end
      end
      return aligned
    end

    def get_all_conserved
      all_conserved = {}
      @a.each_key do |k|
        if @a[k].uniq.size == 1 and 
           @a[k][0] != "-" and 
           @a[k][0] != "*" and 
           @a[k][0] =~ /\S/
          all_conserved[k] = @a[k]
        end
      end
      return all_conserved
    end

  end ## End of HandleAlignment Hash

  class TranscriptCorresponding

    def initialize( obj_handle_alignment, transcript_alignment )
      @ha = obj_handle_alignment
      @ta = transcript_alignment
    end

    def check_correspondence
      @ha.gids.each_key do |gid|
        roop_size = ( ( @ta[gid].size / 3 ) - 1 ) 
        for i in ( 1..roop_size )
          next unless @ha.cp[i][ @ha.gids[gid] ]
          cp_pos = @ha.cp[i][ @ha.gids[gid] ] - 1
          triplet4t = ""
          triplet4t = @ta[gid][(cp_pos*3)] + @ta[gid][(cp_pos*3+1)] + @ta[gid][(cp_pos*3+2)]
          triplet4p = @ha.c[ (i) ][ @ha.gids[gid] ]
p [triplet4t, triplet4p]
        end
      end
      return true
    end

    def show_corresponding_triplet( gid, aa_pos )
      corresp_triplet = ""
      ts_pos = ( @ha.cp[aa_pos][ @ha.gids[gid] ] * 3 ) - 2
      corresp_triplet = @ta[gid][ts_pos] + @ta[gid][(ts_pos+1)] + @ta[gid][(ts_pos+2)]
      return corresp_triplet
    end

  end

  def CreateAlignmentHash( alignment_file_in_pir_format, pirformat=true )

    ch = {}; h  = {}; defline = "";
    ff = FlatFile.new(FastaFormat, open(alignment_file_in_pir_format) )
    ff.each do |e|
      if pirformat == true
        defline = e.definition.slice(/\S+\;(\S+)/, 1)
      else
        defline = e.definition
      end
      if pirformat == true
        h[ defline ] = e.aaseq.split("")
      else
        h[ defline ] = e.naseq.split("")
      end
    end
    ## Check #########
    if pirformat == true and h.values.collect{|v| v.size }.uniq.size > 1
      print "The size of each sequence in the alignment is not same\n"; exit;
    end
    ################
    return h
 
  end

  def CreateConsensus( alignment_file_in_pir_format )

    ch = {}
    h = CreateAlignmentHash( alignment_file_in_pir_format, true )
    ch[0] = []
    h.each_key{|gid| ch[ 0 ] << gid }
    alignment_length = h.values[0].size
    for i in 1..alignment_length
      ch[ i ] = []
      h.each_key do |defline|
        ch[ i ] << h[ defline ][i]
      end
    end
    return ch

  end

#  module_function  :CreateConsensus
#  module_function  :CreateAlignmentHash

end

if $0 == __FILE__

  sample_aa_file = "./sample_files/factory/all_bra.cand.pir"
  ch = ParseAlign.CreateConsensus( sample_aa_file )
#  ch = ParseAlign.CreateConsensus(ARGV.shift)
  pa = ParseAlign::HandleAlignment.new( ch )

  pa.s.each_key do |num|
    print num, "\t"
    p pa.s[num]
    p ch[ num ]
    p pa.cp[num]
  end

end
