#!/usr/bin/ruby

require 'rubygems'
require 'bio'
include Bio


module ParseAlign

  class HandleAlignment

    def initialize( consensus_hash )
      @c = consensus_hash
      @a = get_aligned
      @s = get_stead
      @cp = get_compensation
    end

    attr_reader :c, :a, :s, :cp

    def get_compensation
      cp = {}
      sorted_c_keys = @c.keys.sort{|x, y| x.to_i <=> y.to_i }
      sorted_c_keys.each{|k| cp[k] = [] }
      num_of_species = @c.values[0].size
      non_gap_pos = []
      for j in 0..( num_of_species - 1 )
        non_gap_pos = sorted_c_keys.collect{|k| k if @c[k][j] != "-" and @c[k][j] != "*" }.compact.sort
        non_gap_pos.each_with_index do |k, m|
          cp[k][j] = m
        end
      end
      return cp
    end

    def get_stead
      stead = {}
      @a.each_key do |k|
        if @a[k].uniq.size == 1 and 
           @a[k][0] != "-" and 
           @a[k][0] != "*" and 
           @a[k][0] =~ /\S/
          stead[k] = @a[k]
        end
      end
      return stead
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
      @c.each_key do |k|
      end
    end

  end ## End of Create Alignment Hash

  def CreateAlignment( alignment_file_in_pir_format )

    ch = {}; h  = {}; defline = "";
    ff = FlatFile.new(FastaFormat, open(alignment_file_in_pir_format) )
    ff.each do |e|
      defline = e.definition.slice(/\S+\;(\S+)/, 1)
      h[ defline ] = e.naseq.split("")
    end
    ## Check #########
    if h.values.collect{|v| v.size }.uniq.size > 1
      print "The size of each sequence in the alignment is not same\n"; exit;
    end
    ################
    return h
 
  end

  def CreateConsensus( alignment_file_in_pir_format )

    ch = {}
    h = CreateAlignment( alignment_file_in_pir_format )
    alignment_size = h.values[0].size
    for i in 0..(alignment_size - 1)
      ch[ i ] = []
      h.each_key do |defline|
        ch[ i ] << h[ defline ][i]
      end
    end
    return ch

  end

  module_function  :CreateConsensus
  module_function  :CreateAlignment

end

if $0 == __FILE__

  ch = ParseAlign.CreateConsensus(ARGV.shift)
  pa = ParseAlign::HandleAlignment.new( ch )

  pa.s.each_key do |num|
    print num, "\t"
    p pa.s[num]
    p ch[ num ]
    p pa.cp[num]
  end


end
