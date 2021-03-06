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
    sample_ts_file = "./sample_files/factory/all_bra.cand.nt.fa"
    @obj = ClassParseAlign.new
    @cah = @obj.CreateAlignmentHash(sample_aa_file)
    @cth = @obj.CreateAlignmentHash(sample_ts_file, false )
    @ch  = @obj.CreateConsensus( sample_aa_file )
    @ha  = ClassParseAlign::HandleAlignment.new( @ch )
    @tc  = ClassParseAlign::TranscriptCorresponding.new( @ha, @cth  )
  end

  def test_foo
    assert_equal("foo", @obj.huga  )
  end

  def test_CreateAlignmentHash
    assert_equal( Hash, @cah.class  )
    assert_equal( 11, @cah.keys.size  )
    assert_equal( Hash, @cth.class  )
  end

  def test_get_compensation
    assert_equal( Hash, @ha.get_compensation.class )
    assert_equal( nil, @ha.get_compensation[1][0] )
    assert_equal( nil, @ha.get_compensation[1][0] )
    assert_equal( 1, @ha.get_compensation[1][8] )
    assert_equal( nil, @ha.get_compensation[563][0] )
    assert_equal( 482, @ha.get_compensation[564][10] )
  end

  def test_get_aligned
    assert_equal( Hash, @ha.get_aligned.class )
    assert_equal( 11, @ha.get_aligned.values[0].size )
  end

  def test_get_all_conserved
    assert_equal( Hash, @ha.get_all_conserved.class )
    assert_equal( 11, @ha.get_all_conserved.values[0].size )
    assert_equal( ["M"], @ha.get_all_conserved.values[0].uniq )
    assert_equal( ["L"], @ha.get_all_conserved[366].uniq )
  end

  def test_check_correspondence
    assert_equal( true, @tc.check_correspondence )
  end

  def test_show_corresponding_triplet_by_protein_position
    gid = "ENSMUSP00000027859"
    pos_aa = 90
    assert_equal(Bio::Sequence::NA, @tc.show_corresponding_triplet_by_protein_position(gid, pos_aa).class)
    assert_equal("atg", @tc.show_corresponding_triplet_by_protein_position(gid, pos_aa))
    pos_aa = 91
    assert_equal("ata", @tc.show_corresponding_triplet_by_protein_position(gid, pos_aa))
    pos_aa = 92
    assert_equal("gtg", @tc.show_corresponding_triplet_by_protein_position(gid, pos_aa))
    pos_aa = 93
    assert_equal("aca", @tc.show_corresponding_triplet_by_protein_position(gid, pos_aa))
  end

end



