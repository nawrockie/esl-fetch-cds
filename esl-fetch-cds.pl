#!/usr/bin/env perl
# 
# esl-fetch-cds.pl: fetch a CDS given a NCBI-format specification of 
#                   the exons.
# EPN, Tue Mar 24 09:02:48 2015
# 
# This script uses BioEasel's SqFile module and will create many
# temporary fasta files and .ssi indexes of them, deleting each
# one after it is done using it. 

use strict;
use Getopt::Long;
use Bio::Easel::SqFile;

my $in_cfile     = "";    # input name of input file to split up, 1st cmd line arg
my $outfile_root = undef; # root for name of output file, default is $in_cfile, changed if -oroot used
my $outfile_dir  = "";    # dir for output files, pwd unless -odir is used   
my $idfetch      = "/netopt/ncbi_tools64/bin/idfetch";
my $fasta_linelen = 80; # TODO: make this settable at command line
my $outdir       = undef;

&GetOptions( "odir=s" => \$outdir );

my $usage;
$usage  = "esl-fetch-cds.pl [OPTIONS] <input coordinate file>\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t-odir <s>: write temporary files to directory <s>, instead of cwd\n";

if(scalar(@ARGV) != 1) { die $usage; }
($in_cfile) = @ARGV;

# validate input args
if(! -e $in_cfile) { die "ERROR $in_cfile does not exist"; }
if(defined $outdir) { 
  if(! -d $outdir) { die "ERROR, with -odir <s>, directory <s> must exist, it does not"; }
  $outdir =~ s/\/*//; # remove trailing slash if it's there
  $outdir .= "/"; # add trailing slash
}

my @orig_accn_A   = (); # array of original (protein) accessions for each CDS, which we'll use to name output CDSs
my @codon_start_A = (); # array of 'codon_start' for each CDS, if absent for a CDS, default value is '1'
my @cds_info_AA   = (); # 2D array: each element is an arrays that specifies the information for a CDS, 
                        # and is an array itself of elements of this format: <accession.version>:<start>:<stop>:<strand>
                        # each of which describes an exon.

# parse the coordinate file
parseCoordsFile($in_cfile, \@orig_accn_A, \@codon_start_A, \@cds_info_AA);

#debugPrintCdsInfo(\@orig_accn_A, \@cds_info_AA);

# for each CDS, fetch all sequences for that CDS using idfetch
# POSSIBLE TODO: use idfetch to get sequence length (if possible) then fetch as many seqs as possible
# until we get to a certain max length, then fetch all current CDSs
# TODO: reorganize this code, put in one or two subroutines.
my $ncds = scalar(@orig_accn_A);
for(my $i = 0; $i < $ncds; $i++) { 
  my $orig_accn = $orig_accn_A[$i];
  my $codon_start  = $codon_start_A[$i];
  my $tmp_seqfile  = "tmp.$orig_accn.fa";
  my $tmp_accnfile = "tmp.$orig_accn.accn";
  if(defined $outdir) { 
    $tmp_accnfile = $outdir . $tmp_accnfile;
    $tmp_seqfile  = $outdir . $tmp_seqfile;
  }
  my @unlink_A = (); # files to unlink
  push(@unlink_A, ($tmp_seqfile, $tmp_accnfile));
  my $nexon = scalar(@{$cds_info_AA[$i]});
  my $nseq  = 0; 
  my @seq_A = ();
  my @fetch_info_AA = (); # information for Bio::Easel::SqFile->fetch_subseqs()
                          # elements are arrays with 4 elements: 
                          # <new-name>, <start>, <stop>, <source-name>

  # for each segment/exon
  for(my $e = 0; $e < $nexon; $e++) { 
    my ($seq, $start, $stop, $strand) = split(":", $cds_info_AA[$i][$e]);
    # check if we need to fetch this sequence 
    my $do_fetch = 1; # yes, until proven otherwise
    for (my $s = 0; $s < scalar(@seq_A); $s++) { 
      if($seq_A[$s] eq $seq) { 
        $do_fetch = 0; 
        $s = scalar(@seq_A); # breaks innermost for loop
      }
    }
    push(@seq_A, $seq);
    if($do_fetch) { # we need to fetch this sequence
      # first, we use idfetch to get the seqid for this accession
      my $output_char = (scalar(@seq_A) == 1) ? ">" : ">>";
      my $cmd = "echo $seq > $tmp_accnfile";
      runCommand($cmd);
      # as we do the idfetch, rename the sequence to $seq
      my $cmd = "$idfetch -t 5 -c 1 -G $tmp_accnfile | sed 's/^>\\S*/>$seq/' $output_char $tmp_seqfile";
      runCommand($cmd);
      # TODO: have a command line option for appending commands to a log file, and update runCommand() to output to that file
    }
  }
  # we have all source sequences for all of the exons in $tmp_seqfile, 
  # index it and fetch the subsequences we want
  my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $tmp_seqfile });
  push(@unlink_A, $tmp_seqfile . ".ssi"); # it's very impt to remove this when we're done with it
  my @fetch_info_AA = (); # information for Bio::Easel::SqFile->fetch_subseqs()
                          # elements are arrays with 4 elements: 
                          # <new-name>, <start>, <stop>, <source-name>
  my $cds_name = $orig_accn . ":codon_start" . $codon_start . ":";
  my $expected_strand = ""; # strand of first exon
  my $expected_seq    = ""; # sequence name of first exon
  my $have_multiple_seqs = 0; # set to '1' if we are fetching from more than one source sequence
  for(my $e = 0; $e < $nexon; $e++) { 
    my ($seq, $start, $stop, $strand, $ic_start, $ic_stop) = split(":", $cds_info_AA[$i][$e]);
    my $start2print = ($ic_start) ? "<" . $start : $start;
    my $stop2print  = ($ic_stop)  ? ">" . $stop  : $stop;
    # populate @fetch_info_AA with necessary information
    my $newname = $orig_accn . ":" . $start2print . ":" . $stop2print . ":" . $strand;
    $cds_name  .= $seq . ":" . $start2print . ":" . $stop2print . ":" . $strand . ":";

    # make sure that we have the same strand for all exons, 
    # and keep track of whether we have multiple seqs
    if($e == 0) { 
      $expected_strand = $strand; 
      $expected_seq    = $seq;
    }
    else { 
      if($strand ne $expected_strand) {
        die "ERROR unexpectedly got different strands for different exons of same exon (line: $i+1)";
      }
      if($seq ne $expected_seq) { 
        $have_multiple_seqs = 1;
      }
    }

    if($strand eq "+") { 
      push(@fetch_info_AA, [$newname, $start, $stop, $seq]);
    }
    else { 
      push(@fetch_info_AA, [$newname, $stop, $start, $seq]);
    }
  }

  if($expected_strand eq "-") { # reverse strand, reverse the exons
    if($have_multiple_seqs) { 
      # ERROR, we want to fetch from multiple source sequences on the opposite strand,
      # there's no way to know what the proper order of exons is, is it reversed
      # like it would be with a single source sequence? Or not. If we ever encounter this 
      # then we'll have to figure it out.
      die "ERROR unexpectedly have multiple source sequences with exons on the negative strand. The code 'thought' this was impossible.";
    }
    @fetch_info_AA = reverse @fetch_info_AA;
  }

  my $cds_fa = $sqfile->fetch_subseqs(\@fetch_info_AA, -1, undef); # -1 makes it infinite length
  # we now have the CDS sequence in memory, remove all the header lines,
  $cds_fa =~ s/\n*\>.+\n//g;
  # now output it in chunks of $fasta_linelen
  outputSequence(">" . $cds_name, $cds_fa, $fasta_linelen, undef);

  foreach my $file (@unlink_A) { 
    if(-e $file) { 
      unlink $file;
    }
  }
  $sqfile->close_sqfile();
} # end of 'for(my $i = 0; $i < $ncds; $i++)'    

exit 0;

##############
# SUBROUTINES 
##############
#
# Subroutine: parseCoordsFile()
# Args:       $cfile:          path to coordinate file
#             $orig_accn_AR:   filled here; ref to array of old accessions, one per CDS
#             $codon_start_AR: filled here; ref to array of codon_start, one per CDS
#             $cds_info_AAR:   filled here; 2D array of information about each CDS to fetch
#                              each element includes the information for one CDS and
#                              is itself an array of strings in the following format:
#                              <accession.version>.<start>.<stop>.<strand>
#                              where <start> <= <stop> and <strand> is either '+' or '-'.
#                            
# Returns:    void
# Dies:       if unable to parse the file

sub parseCoordsFile {
  my $sub_name = "ParseCoordsFile()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cfile, $orig_accn_AR, $codon_start_AR, $cds_info_AAR) = (@_);

  my ($orig_line, $line, $strand, $orig_accn); 

  open(IN, $cfile) || die "ERROR unable to open file $cfile for reading";
  while($line = <IN>) { 
    my $orig_line = $line; # keep this saved so we can output it in case of an ERROR
    my @cur_cds_A = ();
    chomp $line;

    # Determine if we have 'codon_start' information, or not, and remove it if we do
    my $codon_start = 1; # default value
    if($line =~ s/\t([123])$//) { 
      $codon_start = $1;
    }
    
    # Determine if 'complement' exists (which gives us strand information) and get 
    # rid of it if it does.
    # examples: 
    # AJJ35879.1	complement(CP009997.1:2173412..2176090)
    # EEH43584.1	complement(join(KN275973.1:226623..226774, KN275973.1:226854..229725))
    # AFJ92633.1	JQ690867.1:<1..400	2
    # BAJ21116.1	AB590961.1:10..>891
    # AAV98535.1	AY680453.1:<1..>489
    
    if($line =~ s/complement\(//) { # negative strand
      $strand = "-";
      if($line =~ s/\)$//) { ; } # remove final character, it should be a ')' matching the '(' in 'complement(' we just removed
      else { die "ERROR unable to parse (failure to find matching ) for 'complement ('; line: $orig_line"; }
    }
    else { # positive strand
      $strand = "+";
    }
    
    # Next deal with lines with 'join' indicating more than one subsequence (i.e. multiple exons)
    # example: 
    # EEH43584.1	join(KN275973.1:226623..226774, KN275973.1:226854..229725)
    
    my ($accn, $start, $stop, $ic_start, $ic_stop);
    if($line =~ s/^(\S+)\s+join\(//) { # get rid of everything up to 'join('
      $orig_accn = $1;
      # now: KN275973.1:226623..226774, KN275973.1:226854..229725))
      chomp $line;
      if($line =~ s/\)$//) { # remove trailing '))'
        # now: KN275973.1:226623..226774, KN275973.1:226854..229725
        # process each /,\s*/ seperated element at a time
        my @elA = split(/\s*\,\s*/, $line);
        my $nel = scalar(@elA);
        my $start_stop = ""; # we'll append each exon to this 
        my $nexon = 0;
        for(my $i = 0; $i < $nel; $i++) { 
          ($accn, $start, $stop, $ic_start, $ic_stop) = tokenBreakdown($elA[$i], $orig_line);
          push(@cur_cds_A, ("$accn:$start:$stop:$strand:$ic_start:$ic_stop"));
        }
      } # end of 'if($line =~ s/\)$//)' 
      else { 
        die "ERROR unable to parse (failure to find matching ) for 'join ('; line: $orig_line";
      }
    } # end of 'if($line =~ s/^\S+\s+\(join\(//) {'
    elsif($line =~ /^(\S+)\s+.+\:\<?\d+\.\.\>?\d+$/) { # no 'join' just one subsequence, easy case
      $orig_accn = $1;
      #AIV71043.1	CP009235.1:2698208..2701135
      $line =~ s/^\S+\s+//; 
      #now: CP009235.1:2698208..2701135
      ($accn, $start, $stop, $ic_start, $ic_stop) = tokenBreakdown($line, $orig_line);
      push(@cur_cds_A, ("$accn:$start:$stop:$strand:$ic_start:$ic_stop"));
    }
    else { 
      die "ERROR unable to parse line: $orig_line";
    }
    # done parsing this line
    push(@{$orig_accn_AR}, $orig_accn); 
    push(@{$codon_start_AR}, $codon_start);
    push(@{$cds_info_AAR}, [@cur_cds_A]);
  } # end of 'while($line = <IN>)'
  return
}  

# Subroutine: tokenBreakdown()
# Args:       $token:        path to coordinate file
#             $orig_line:    original line, only used if there's an error to inform
#                            user about which line the error occurred on.
# Return:      5 values:
#              $accn:       the accession
#              $start:      the start coordinate
#              $stop:       the stop coordinate
#              $ic_start:   '1' if a '<' character occurs before the
#                           start coordinte indicating an incomplete
#                           CDS at the 'start' coordinate (note this
#                           could be either the 5' or 3' end (!)
#                           depending on the strand.
#              $ic_stop:    '1' if a '>' character occurs before the stop
#                           coordinate indicating an incomplete CDS at
#                           the 'stop' coordinate (note this could be
#                           either the 5' or 3' end (!) depending on
#                           the strand.
#
# Dies:       if unable to parse the token
sub tokenBreakdown { 
#
# breakdown a token like this:
#
# CP009235.1:2698208..2701135
# JGZA01000019.1:<1..1949 # flush with beginning of sequence
# GQ358600.1:<1..>453     # flush with end of sequence
# 
# into just <accn> <start> <stop>
# and determine if the subsequence ends (start and/or stop) are
# flush with the beginning or end of the source sequence.
# as indicated by a '<' before the start coordinate and 
# a '>' before the end coordinate.
#
# Return values for examples above:
# 
# CP009235.1 2698208 2701135 0 0 
# JGZA01000019.1 1 1949 1 0 
# GQ358600.1:1 453 1 1
# 
  if(scalar(@_) != 2) { die "ERROR entered tokenBreakdown() with wrong number of args"; }
  my ($token, $orig_line) = (@_);
  my ($ic_start, $ic_stop);
  
  if($token =~ /^(.+)\:(\<?\d+)\.\.(\>?\d+)/) { 
    my ($accn, $start, $stop) = ($1, $2, $3);
    if($start =~ s/^\<//) { $ic_start = 1; }
    else                  { $ic_start = 0; }
    if($stop  =~ s/^\>//) { $ic_stop  = 1; }
    else                  { $ic_stop  = 0; }
    if($start > $stop) { die "ERROR unexpectedly start > stop in $token on line: $orig_line"; }
    return ($accn, $start, $stop, $ic_start, $ic_stop);
  }
  else { 
    die "ERROR unable to breakdown $token, in line: $orig_line"; 
  }
}

# Subroutine: debugPrintCdsInfo()
# Args:       $orig_accn_AR: ref to array of old accessions, one per CDS
#             $cds_info_AAR: 2D array of information about each CDS to fetch
# Return:     void
sub debugPrintCdsInfo { 
  if(scalar(@_) != 2) { die "ERROR entered tokenBreakdown() with wrong number of args"; }
  my ($orig_accn_AR, $cds_info_AAR) = (@_);

  my $ncds = scalar(@{$orig_accn_AR});
  for(my $i = 0; $i < $ncds; $i++) { 
    printf("$orig_accn_AR->[$i]");
    my $nexon = scalar(@{$cds_info_AAR->[$i]});
    for(my $e = 0; $e < $nexon; $e++) { 
      printf(" $cds_info_AAR->[$i][$e]");
    }
    printf("\n");
  }
  return 0;
}

# Subroutine: runCommand()
# Args:       $cmd:            command to run, with a "system" command;
#
# Returns:    void
# Dies:       if $cmd fails

sub runCommand {
  my $sub_name = "runCommand()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cmd) = @_;

  system($cmd);

  if($? != 0) { die "ERROR command failed:\n$cmd\n"; }

  return;
}

# Subroutine: outputSequence()
# Args:       $header:  fasta header line
#             $seq:     actual sequence
#             $linelen: length for each line
#             $FH:      file handle to print to, if undef print to stdout

# Returns:    void

sub outputSequence {
  my $sub_name = "outputSequence()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($header, $seq, $linelen, $FH) = @_;
  chomp $seq;

  if(defined $FH) { print $FH $header . "\n"; }
  else            { print $header . "\n";     }
  
  my $length = length($seq);
  my $p = 0;
  while($p < $length) { 
    my $toprint = substr($seq, $p, $linelen);
    if(defined $FH) { print $FH $toprint . "\n"; }
    else            { print $toprint . "\n";     }
    $p += $linelen;
  }

  return;
}
