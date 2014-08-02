#!/usr/bin/ruby

require 'test/unit'
require './GFF3Parser'

class ClassGFF3Parser
  include GFF3Parser
end

class Test_GFF3Parser < Test::Unit::TestCase

  def setup
      @obj = ClassGFF3Parser.new
  end

  def test_foo
    assert_equal("foo", @obj.create_gene_hash )
  end

end



