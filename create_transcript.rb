

require 'rubygems'
require 'bio'
include Bio

def create_revtrans( aas, defval )
  revtrans_seq = ""
  table = CodonTable[1]
  aas.split("").each do |aa|
    if table.revtrans( aa ) == nil
      revtrans_seq += "nnn"
    else
      revtrans_seq += table.revtrans( aa ).sample
    end
  end
  return Sequence::NA.new( revtrans_seq )
end

ff = FlatFile.new( FastaFormat, ARGF )
ff.each do |e|
#   print e.aaseq.to_fasta( e.definition, 60 )
  print create_revtrans( e.aaseq, e.definition ).to_fasta( e.definition, 60 )
end


