#!/usr/bin/ruby

require 'rubygems'
require 'minitest/autorun'
require './ParseAlign'

MiniTest.autorun

class ClassParseAlign
  include ParseAlign
end

class Test_ParseAlign < MiniTest::Test

  def setup
    sample_aa_file = "./sample_files/factory/all_bra.cand.pir"
    @obj = ClassParseAlign.new
    @cah = @obj.CreateAlignmentHash(sample_aa_file)
    @ch  = @obj.CreateConsensus( sample_aa_file )
    @ha  = ClassParseAlign::HandleAlignment.new( @ch )
  end

  def test_foo
    assert_equal("foo", @obj.huga  )
  end

  def test_CreateAlignmentHash
    assert_equal( Hash, @cah.class  )
    assert_equal( 11, @cah.keys.size  )
  end

  def test_get_compensation
    assert_equal( Hash, @ha.get_compensation.class )
    assert_equal( nil, @ha.get_compensation[0][0] )
    assert_equal( nil, @ha.get_compensation[0][0] )
    assert_equal( 0, @ha.get_compensation[0][8] )
    assert_equal( nil, @ha.get_compensation[562][0] )
    assert_equal( 480, @ha.get_compensation[562][10] )
  end

  def test_get_aligned
    assert_equal( Hash, @ha.get_aligned.class )
    assert_equal( 11, @ha.get_aligned.values[0].size )
  end

  def test_get_all_conserved
    assert_equal( Hash, @ha.get_all_conserved.class )
    assert_equal( 11, @ha.get_all_conserved.values[0].size )
    assert_equal( ["m"], @ha.get_all_conserved.values[0].uniq )
    assert_equal( ["l"], @ha.get_all_conserved[366].uniq )
  end

end

#  def CreateAlignmentHash( alignment_file_in_pir_format, pirformat="true" )


