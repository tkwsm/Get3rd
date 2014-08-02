#!/usr/bin/ruby

require 'rubygems'
require 'bio'
include Bio

class ParseAlign

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
    return aligned
  end

  def get_all_conserved
    @c.each_key do |k|
    end
  end

end

ch = {}
h  = {}
ff = FlatFile.new(FastaFormat, ARGF)
defline = ""
ff.each do |e|
  defline = e.definition.slice(/\S+\;(\S+)/, 1)
  h[ defline ] = e.naseq.split("")
end

### Check #########
  if h.values.collect{|v| v.size }.uniq.size > 1
    print "The size of each sequence in the alignment is not same\n"; exit;
  end
#################

alignment_size = h.values[0].size
for i in 0..(alignment_size - 1)
  ch[ i ] = []
  h.each_key do |defline|
    ch[ i ] << h[ defline ][i]
  end
end

pa = ParseAlign.new( ch )



