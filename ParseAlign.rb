#!/usr/bin/ruby

require 'rubygems'
require 'bio'
include Bio


module ParseAlign

  class HandleAlignment

    def initialize( consensus_hash )
      @c = consensus_hash
      @a = get_aligned
    end

    attr_reader :c, :a

    def get_aligned
      aligned = {}
      @c.each_key do |k|
        if @c[k].include?("-")
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

pa = ParseAlign::HandleAlignment.new(ParseAlign.CreateConsensus(ARGV.shift)).c
pa.each_key do |num|
  print num, "\t"
  p pa[num]
end



