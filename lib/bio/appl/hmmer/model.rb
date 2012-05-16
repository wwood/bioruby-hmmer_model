## Example HMM model file, taken from PFAM 26.

# HMMER3/b [3.0 | March 2010]
# NAME  1-cysPrx_C
# ACC   PF10417.4
# DESC  C-terminal domain of 1-Cys peroxiredoxin
# LENG  40
# ALPH  amino
# RF    no
# CS    yes
# MAP   yes
# DATE  Mon Sep 26 01:36:52 2011
# NSEQ  86
# EFFN  26.221497
# CKSUM 2893062708
# GA    20.40 20.40;
# TC    20.40 20.50;
# NC    20.30 20.30;
# BM    hmmbuild HMM.ann SEED.ann
# SM    hmmsearch -Z 15929002 -E 1000 --cpu 4 HMM pfamseq
# STATS LOCAL MSV       -7.4458  0.71948
# STATS LOCAL VITERBI   -7.8857  0.71948
# STATS LOCAL FORWARD   -4.3017  0.71948
# HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
#             m->m     m->i     m->d     i->m     i->i     d->m     d->d
#   COMPO   2.25274  4.33630  2.74834  2.65826  3.87771  2.70273  3.95751  3.25125  2.56848  2.82101  4.06536  3.21194  2.50202  2.97228  3.39798  2.99665  2.70159  2.62185  3.52465  3.78187
#           2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
#           0.00153  6.88065  7.60300  0.61958  0.77255  0.00000        *
#       1   0.35002  6.53529  7.25102  7.30399  7.52698  2.51105  7.66616  7.05965  7.14684  6.75377  3.79371  6.32251  6.23186  7.14685  7.02553  1.76469  5.22656  6.02391  8.86510  7.90937      1 - H
#           2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
#           0.00153  6.88065  7.60300  0.61958  0.77255  0.48576  0.95510
#       2   3.84702  5.95842  6.58626  5.97181  2.06407  5.79758  6.11757  2.43011  5.76033  0.58540  3.17079  5.96610  6.11889  5.84497  5.74893  5.11667  4.83456  2.97676  3.47925  3.13829      2 - H
#           2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
#           0.00153  6.88065  7.60300  0.61958  0.77255  0.48576  0.95510

module Bio
  class Hmmer
    class Model
      class ParseException<Exception; end

      TAGS = %w(NAME ACC DESC LENG ALPH RF CS MAP DATE NSEQ EFFN CKSUM GA TC NC BM SM STATS_LOCAL_MSV)
      TAGS.each do |info|
        attr_accessor info.downcase.to_sym
      end

      # An array of 4 element arrays from STATS lines (e.g. STATS LOCAL FORWARD   -4.3017  0.71948) => ['LOCAL','FORWARD',-4.3017,0.71948]
      attr_accessor :stats

      # Letters used in the model (amino acids). The order of this array is the same as the order
      attr_accessor :alphabet

      # ["m->m", "m->i", "m->d", "i->m", "i->i", "d->m", "d->d"]
      attr_accessor :transitions

      # Probabilities for matches in the model (first line of each position in the main model part of the HMMER3/b file format)
      # You probably want to use #match_probability not read this variable directly
      attr_accessor :match_emissions

      # Probabilities for inserts in the model (second line of each position in the main model part of the HMMER3/b file format)
      attr_accessor :insert_emissions

      # Probabilities for state transitions in the model (third line of each position in the main model part of the HMMER3/b file format)
      attr_accessor :state_transitions

      # Return the probability given in the HMM of the match probability of a particular position.
      #
      # Note that the hmm_position is 1-based - the first position is 1, not 0
      #
      # Assumes that the letter supplied is in the HMM's alphabet
      def match_probability(hmm_position, letter)
        index = @alphabet.index(letter)
        match_emissions[hmm_position+1][index]
      end

      # Parse a HMM file fed in through some IO, and return an instantiated
      # Bio::HMMER::Model with all the informations filled in
      #
      # e.g. PF10417.4 has length 40 (as of PFAM v26):
      #
      #    Bio::Hmmer::Model.parse(File.open('test/data/PF10417.4.hmm')).leng => 40
      def self.parse(io)
        model = Bio::Hmmer::Model.new

        lines = io.readlines

        unless lines[0].match(/^HMMER3\/b/)
          raise ParseException, "Only HMMER3/b files are currently supported, sorry. The first line of the HMMER model file causing the problem is #{lines[0]}"
        end
        lineno = 1

        # Majority of tags at the top of the model
        while TAGS.include?(lines[lineno].split(/ +/)[0])
          splits = lines[lineno].strip.split(/ +/)
          lit = splits[0].downcase
          tag = TAGS.select{|t| t==lit}[0]
          value = splits[1..(splits.length-1)].join(' ')

          # cast appropriately
          value = value.to_i if %w(leng nseq cksum).include?(lit)
          value = value.to_f if %w(effn).include?(lit)
          value = [splits[1].to_f, splits[2].to_f] if %w(ga tc nc).include?(lit)

          model.send("#{lit}=".to_sym,value)
          lineno += 1
        end

        # STATS_LOCAL_MSV STATS_LOCAL_VITERBI STATS_LOCAL_FORWARD
        while matches = lines[lineno].match(/^STATS ([^ ]+) ([^ ]+) +([\-\d\.]+) +([\-\d\.]+)$/)
          model.stats ||= []
          model.stats.push [matches[1], matches[2], matches[3].to_f, matches[4].to_f]
          lineno += 1
        end

        # Next expect HMM   A C D ...
        splits = lines[lineno].strip.split(/ +/)
        unless splits[0]=='HMM'
          raise ParseException, "Unexpected HMM file line encountered, expected the beginning of the main model section e.g. HMM          A        C        D        E ...: #{lines[lineno]}"
        end
        model.alphabet = splits[1..(splits.length-1)]
        lineno += 1

        # m->m     m->i     m->d     i->m     i->i     d->m     d->d
        splits = lines[lineno].strip.split(/ +/)
        unless splits.length == 7
          raise ParseException, "Unexpected line in HMM file, expected '            m->m     m->i     m->d     i->m     i->i     d->m     d->d', found #{lines[lineno]}"
        end
        model.transitions = splits
        lineno += 1

        #skip the next lines since my application does not require these
        splits = lines[lineno].strip.split(/ +/)
        lineno += 1 if splits[0]=='COMPO'
        lineno += 2

        model.match_emissions = []
        model.insert_emissions = []
        model.state_transitions = []

        # Iterate through the model, one line at a time
        model_iterate = 0
        while !(lines[lineno].match(/^\/\//))
          splits = lines[lineno].strip.split(/ +/)

          unless splits.length > model.alphabet.length
            raise ParseException, "Unexpected line format encountered in the main model: #{lines[lineno]}"
          end

          # Check that the state numbers are in order
          model_iterate += 1
          unless splits[0].to_i == model_iterate
            raise ParseException, "Unexpected model state number: #{splits[0]}"
          end
          probabilities = splits[1..model.alphabet.length].collect do |s|
          # From the HMMER manual:
          # All probability parameters are all stored as negative natural log probabilities with
          # five digits of precision to the right of the decimal point, rounded. For example,
          # a probability of 0.25 is stored as − log 0.25 = 1.38629.
            10.0**-(s.to_f)
          end
          model.match_emissions.push probabilities
          lineno += 1

          splits = lines[lineno].strip.split(/ +/)
          unless splits.length == model.alphabet.length
            raise ParseException, "Unexpected line format encountered in the main model: #{lines[lineno]}"
          end
          model.insert_emissions.push splits.collect{|s| 10.0**-(s.to_f)}
          lineno += 1

          splits = lines[lineno].strip.split(/ +/)
          unless splits.length == model.transitions.length
            raise ParseException, "Unexpected line format encountered in the main model: #{lines[lineno]}"
          end
          model.state_transitions.push splits.collect{|s| 10.0**-(s.to_f)}
          lineno += 1
        end

        return model
      end
    end
  end
end