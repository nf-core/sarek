=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>
    
=cut

=head1 NAME

 GeneSplicer

=head1 SYNOPSIS

 mv GeneSplicer.pm ~/.vep/Plugins
 ./vep -i variants.vcf --plugin GeneSplicer,[path_to_genesplicer_bin],[path_to_training_dir],[option1=value],[option2=value]

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 runs GeneSplicer (https://ccb.jhu.edu/software/genesplicer/) to get
 splice site predictions.

 It evaluates a tract of sequence either side of and including the
 variant, both in reference and alternate states. The amount of
 sequence included either side defaults to 100bp, but can be modified
 by passing e.g. "context=50" as a parameter to the plugin.

 Any predicted splicing regions that overlap the variant are reported
 in the output with one of four states: no_change, diff, gain, loss

 There follows a "/"-separated string consisting of the following data:

 1) type (donor, acceptor)
 2) coordinates (start-end)
 3) confidence (Low, Medium, High)
 4) score

 Example: loss/acceptor/727006-727007/High/16.231924

 If multiple sites are predicted, their reports are separated by ",".

 For diff, the confidence and score for both the reference and alternate
 sequences is reported as REF-ALT.

 Example: diff/donor/621915-621914/Medium-Medium/7.020731-6.988368

 Several parameters can be modified by passing them to the plugin string:

 context    : change the amount of sequence added either side of
              the variant (default: 100bp)
 tmpdir     : change the temporary directory used (default: /tmp)
 cache_size : change how many sequences' scores are cached in memory
              (default: 50)

 Example: --plugin GeneSplicer,$GS/bin/linux/genesplicer,$GS/human,context=200,tmpdir=/mytmp

 On some systems the binaries provided will not execute, but can be compiled from source:

   cd $GS/sources
   make
   cd -
   ./vep [options] --plugin GeneSplicer,$GS/sources/genesplicer,$GS/human

 On Mac OSX the make step is known to fail; the genesplicer.cpp file requires modification:

   cd $GS/sources
   perl -pi -e "s/^main  /int main  /" genesplicer.cpp
   make
 

=cut

package GeneSplicer;

use strict;
use warnings;

use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

our %DEFAULTS = (
  context => 100,
  tmpdir  => '/tmp',
  cache_size => 50,
);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);
  
  # we need sequence, so no offline mode unless we have FASTA
  die("ERROR: cannot function in offline mode without a FASTA file\n") if $self->{config}->{offline} && !$self->{config}->{fasta};

  my $params = $self->params;

  my $bin = shift @$params;
  die("ERROR: genesplicer binary not specified\n") unless $bin;
  die("ERROR: genesplicer binary not found\n") unless -e $bin;
  my $test = `$bin 2>&1`;
  die("ERROR: failed to run genesplicer binary:\n$test\n") unless $test =~ /^USAGE/;
  $self->{_bin} = $bin;
  
  my $training_dir = shift @$params;
  die("ERROR: training directory not specified\n") unless $training_dir;
  die("ERROR: training directory not found\n") unless -d $training_dir;
  $self->{_training_dir} = $training_dir;

  # defaults
  $self->{'_param_'.$_} = $DEFAULTS{$_} for keys %DEFAULTS;

  # REST API passes 1 as first param
  shift @$params if $params->[0] && $params->[0] eq '1';

  # set/override with user params
  foreach my $param(@$params) {
    my ($key, $val) = split('=', $param);
    die("ERROR: Failed to parse parameter $param\n") unless defined($key) && defined($val);

    $self->{'_param_'.$key} = $val;
  }

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {
    GeneSplicer => "GeneSplicer predictions"
  };
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;

  # get up and downstream sequences
  my $up_seq = $vf->{slice}->sub_Slice(
    $vf->{start} - $self->{'_param_context'},
    $vf->{start} - 1,
    $vf->strand
  )->seq;

  my $down_seq = $vf->{slice}->sub_Slice(
    $vf->{end} + 1,
    $vf->{end} + $self->{'_param_context'},
    $vf->strand
  )->seq;

  # create ref seq by grabbing reference TVA
  my $ref_seq = join("",
    $up_seq,
    $tva->transcript_variation->get_reference_TranscriptVariationAllele->variation_feature_seq,
    $down_seq
  );

  return {} unless $ref_seq =~ /^[ACGT]+$/;

  # create alt seq
  my $alt_allele = $tva->variation_feature_seq;
  $alt_allele =~ s/\-//g;
  my $alt_seq = $up_seq.$alt_allele.$down_seq;


  return {} unless $alt_seq =~ /^[ACGT]+$/;

  # reverse comp if strands differ
  if($tva->transcript->strand != $vf->strand) {
    reverse_comp(\$ref_seq);
    reverse_comp(\$alt_seq);
  }

  # get results
  my $ref_results = $self->results_from_cache($ref_seq) || $self->results_from_seq($ref_seq);
  my $alt_results = $self->results_from_cache($alt_seq) || $self->results_from_seq($alt_seq);

  # compare results both ways
  my $diff_ref_to_alt = $self->compare_results($ref_results, $alt_results);
  my $diff_alt_to_ref = $self->compare_results($alt_results, $ref_results);

  # get VF pos relative to tested sequence
  my ($vf_start, $vf_end) = ($self->{'_param_context'} + 1, $self->{'_param_context'} + (($vf->{end} - $vf->{start}) + 1));

  # get overlapping losses and gains
  # and map to chromosome coords
  my @losses =
    map {$_->{gl} = 'loss'; $_}
    @{$diff_ref_to_alt->{lost}};

  my @gains = 
    map {$_->{gl} = 'gain'; $_}
    @{$diff_alt_to_ref->{lost}};

  my @diffs = 
    map {$_->{gl} = 'diff'; $_}
    @{$diff_ref_to_alt->{diff}};

  my $return = join(',',
    map {
      join('/',
        $_->[0]->{gl},
        $_->[0]->{type},
        $_->[1]->{end5}.'-'.$_->[1]->{end3},
        $_->[0]->{confidence},
        $_->[0]->{score}
      )
    }
    map {[$_, $self->map_ss_coords($_, $vf)]}
    grep {overlap($vf_start, $vf_end, $_->{end5}, $_->{end3})}
    (@losses, @gains, @diffs)
  );

  # probably of interest to report splice sites were found
  # but no difference between ref and alt
  if(!$return && grep {overlap($vf_start, $vf_end, $_->{end5}, $_->{end3})} @$ref_results) {
    $return = join(',',
      map {
        join('/',
          'no_change',
          $_->[0]->{type},
          $_->[1]->{end5}.'-'.$_->[1]->{end3},
          $_->[0]->{confidence},
          $_->[0]->{score}
        )
      }
      map {[$_, $self->map_ss_coords($_, $vf)]}
      grep {overlap($vf_start, $vf_end, $_->{end5}, $_->{end3})} @$ref_results
    );
  }

  return $return ? { GeneSplicer => $return } : {};
}

sub results_from_seq {
  my $self = shift;
  my $seq = shift;

  # write seqs to file
  my $seq_file = $self->{'_param_tmpdir'}."/genesplicer_$$.fa";
  open SEQ, ">$seq_file" or die("ERROR: Could not write to temporary sequence file $seq_file\n");
  print SEQ ">SEQ\n$seq\n";
  close SEQ;

  my $result_file = $self->{'_param_tmpdir'}."/genesplicer_$$.results";
  
  my $cmd = sprintf(
    '%s %s %s -f %s',
    $self->{'_bin'},
    $seq_file,
    $self->{'_training_dir'},
    $result_file
  );

  my $output = `$cmd 2>&1`;
  unlink($seq_file);

  return [] unless -e $result_file;

  open RES, $result_file;
  my @results;
  
  while(<RES>) {
    chomp;
    my ($end5, $end3, $score, $confidence, $type) = split;

    push @results, {
      end5 => $end5,
      end3 => $end3,
      score => $score,
      confidence => $confidence,
      type => $type
    };
  }
  close RES;

  unlink($result_file);

  push @{$self->{cache}}, { hex => md5_hex($seq), results => \@results};
  shift @{$self->{cache}} while scalar @{$self->{cache}} > $self->{_param_cache_size};

  return \@results;
}

sub results_from_cache {
  my $self = shift;
  my $seq = shift;

  my ($results) = map {$_->{results}} grep {$_->{hex} eq md5_hex($seq)} @{$self->{cache} || []};

  return $results;
}

sub compare_results {
  my $self = shift;
  my $a = shift;
  my $b = shift;

  my (@diff, @lost);

  foreach my $res_a(@$a) {
    my @match = grep {
      $_->{end5} == $res_a->{end5} &&
      $_->{end3} == $res_a->{end3} &&
      $_->{type} eq $res_a->{type}
    } @$b;

    # result not found in b
    if(!@match) {
      push @lost, $res_a;
    }

    # >1 result found
    elsif(scalar @match > 1) {
      warn("WARNING: Found two matches?\n");
    }

    # 1 match
    elsif($match[0]->{score} != $res_a->{score}) {
      my %diff = %$res_a;
      $diff{score}      .= '-'.$match[0]->{score};
      $diff{confidence} .= '-'.$match[0]->{confidence};
      push @diff, \%diff;
    }
  }

  return { diff => \@diff, lost => \@lost};
}

sub map_ss_coords {
  my $self = shift;
  my $res = shift;
  my $vf = shift;

  my $return = {};

  foreach my $coord(qw(end5 end3)) {
    $return->{$coord} = (($res->{$coord} - $self->{'_param_context'}) + $vf->{start}) - 1;
  }

  return $return;
}

1;

