require 'helper'

class TestBioHmmerModel < Test::Unit::TestCase
  should "read a single file OK" do
    hmm = Bio::Hmmer::Model.parse(File.open('test/data/PF10417.4.hmm'))
    assert_equal 86, hmm.nseq
    assert_equal 40, hmm.leng
    assert_in_delta 0.00000005, hmm.match_emissions[0][3], 1e-8
    assert_equal ['LOCAL','FORWARD',-4.3017,0.71948], hmm.stats[2]
  end
  
  should "give the correct match probability" do
    hmm = Bio::Hmmer::Model.parse(File.open('test/data/PF10417.4.hmm'))
    assert_in_delta 0.000204904, hmm.match_probability(3, 'I'), 1e-4
  end
end
