require 'helper'

class TestBioHmmerModel < Test::Unit::TestCase
  should "read a single file OK" do
    hmm = Bio::HMMER::Model.parse(File.open('test/data/PF10417.4.hmm'))
    assert_equal 86, hmm.nseq
    assert_equal 40, hmm.leng
    assert_in_delta 0.000672849, hmm.match_emissions[0][3], 0.000672849/10000
    assert_equal "hmmsearch -Z 15929002 -E 1000 --cpu 4 HMM pfamseq", hmm.sm
    assert_equal ['LOCAL','FORWARD',-4.3017,0.71948], hmm.stats[2]
  end
  
  should "give the correct match probability" do
    hmm = Bio::HMMER::Model.parse(File.open('test/data/PF10417.4.hmm'))
    assert_in_delta 0.025010739, hmm.match_probability(3, 'I'), 0.025010739/10000
  end
end
